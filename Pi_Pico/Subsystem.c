#include "Subsystem.h"
#include "Subsystem_private.h"

/* Block signals (default storage) */
B_Subsystem_T Subsystem_B;

/* Continuous states */
X_Subsystem_T Subsystem_X;

/* Block states (default storage) */
DW_Subsystem_T Subsystem_DW;

/* External inputs (root inport signals with default storage) */
ExtU_Subsystem_T Subsystem_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_Subsystem_T Subsystem_Y;

/* Real-time model */
static RT_MODEL_Subsystem_T Subsystem_M_;
RT_MODEL_Subsystem_T *const Subsystem_M = &Subsystem_M_;
real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
                     uint32_T maxIndex)
{
  real_T frac;
  real_T yL_0d0;
  uint32_T iLeft;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    uint32_T bpIdx;
    uint32_T iRght;

    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  yL_0d0 = table[iLeft];
  return (table[iLeft + 1U] - yL_0d0) * frac + yL_0d0;
}

/* Simplified version of numjac.cpp, for use with RTW. */
void local_numjac( RTWSolverInfo *si, real_T *y, const real_T *Fty, real_T *fac,
                  real_T *dFdy )
{
  /* constants */
  real_T THRESH = 1e-6;
  real_T EPS = 2.2e-16;                /* utGetEps(); */
  real_T BL = pow(EPS, 0.75);
  real_T BU = pow(EPS, 0.25);
  real_T FACMIN = pow(EPS, 0.78);
  real_T FACMAX = 0.1;
  int_T nx = 3;
  real_T *x = rtsiGetContStates(si);
  real_T del;
  real_T difmax;
  real_T FdelRowmax;
  real_T temp;
  real_T Fdiff;
  real_T maybe;
  real_T xscale;
  real_T fscale;
  real_T *p;
  int_T rowmax;
  int_T i,j;
  if (x != y)
    (void) memcpy(x, y,
                  (uint_T)nx*sizeof(real_T));
  rtsiSetSolverComputingJacobian(si,true);
  for (p = dFdy, j = 0; j < nx; j++, p += nx) {
    /* Select an increment del for a difference approximation to
       column j of dFdy.  The vector fac accounts for experience
       gained in previous calls to numjac. */
    xscale = fabs(x[j]);
    if (xscale < THRESH)
      xscale = THRESH;
    temp = (x[j] + fac[j]*xscale);
    del = temp - y[j];
    while (del == 0.0) {
      if (fac[j] < FACMAX) {
        fac[j] *= 100.0;
        if (fac[j] > FACMAX)
          fac[j] = FACMAX;
        temp = (x[j] + fac[j]*xscale);
        del = temp - x[j];
      } else {
        del = THRESH;                  /* thresh is nonzero */
        break;
      }
    }

    /* Keep del pointing into region. */
    if (Fty[j] >= 0.0)
      del = fabs(del);
    else
      del = -fabs(del);

    /* Form a difference approximation to column j of dFdy. */
    temp = x[j];
    x[j] += del;
    Subsystem_step();
    rtsiSetdX(si,p);
    Subsystem_derivatives();
    x[j] = temp;
    difmax = 0.0;
    rowmax = 0;
    FdelRowmax = p[0];
    temp = 1.0 / del;
    for (i = 0; i < nx; i++) {
      Fdiff = p[i] - Fty[i];
      maybe = fabs(Fdiff);
      if (maybe > difmax) {
        difmax = maybe;
        rowmax = i;
        FdelRowmax = p[i];
      }

      p[i] = temp * Fdiff;
    }

    /* Adjust fac for next call to numjac. */
    if (((FdelRowmax != 0.0) && (Fty[rowmax] != 0.0)) || (difmax == 0.0)) {
      fscale = fabs(FdelRowmax);
      if (fscale < fabs(Fty[rowmax]))
        fscale = fabs(Fty[rowmax]);
      if (difmax <= BL*fscale) {
        /* The difference is small, so increase the increment. */
        fac[j] *= 10.0;
        if (fac[j] > FACMAX)
          fac[j] = FACMAX;
      } else if (difmax > BU*fscale) {
        /* The difference is large, so reduce the increment. */
        fac[j] *= 0.1;
        if (fac[j] < FACMIN)
          fac[j] = FACMIN;
      }
    }
  }

  rtsiSetSolverComputingJacobian(si,false);
}                                      /* end local_numjac */

/*
 * This function updates continuous states using the ODE14X fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static int_T rt_ODE14x_N[4] = { 12, 8, 6, 4 };

  time_T t0 = rtsiGetT(si);
  time_T t1 = t0;
  time_T h = rtsiGetStepSize(si);
  real_T *x1 = rtsiGetContStates(si);
  int_T order = rtsiGetSolverExtrapolationOrder(si);
  int_T numIter = rtsiGetSolverNumberNewtonIterations(si);
  ODE14X_IntgData *id = (ODE14X_IntgData *)rtsiGetSolverData(si);
  real_T *x0 = id->x0;
  real_T *f0 = id->f0;
  real_T *x1start = id->x1start;
  real_T *f1 = id->f1;
  real_T *Delta = id->Delta;
  real_T *E = id->E;
  real_T *fac = id->fac;
  real_T *dfdx = id->DFDX;
  real_T *W = id->W;
  int_T *pivots = id->pivots;
  real_T *xtmp = id->xtmp;
  real_T *ztmp = id->ztmp;
  int_T *N = &(rt_ODE14x_N[0]);
  int_T i,j,k,iter;
  int_T nx = 3;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(x0, x1,
                (uint_T)nx*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */

  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  Subsystem_derivatives();
  local_numjac(si,x0,f0,fac,dfdx );
  for (j = 0; j < order; j++) {
    real_T *p;
    real_T hN = h/N[j];

    /* Get the iteration matrix and solution at t0 */

    /* [L,U] = lu(M - hN*J) */
    (void) memcpy(W, dfdx,
                  (uint_T)nx*nx*sizeof(real_T));
    for (p = W, i = 0; i < nx*nx; i++, p++) {
      *p *= (-hN);
    }

    for (p = W, i = 0; i < nx; i++, p += (nx+1)) {
      *p += 1.0;
    }

    rt_lu_real(W, nx,
               pivots);

    /* First Newton's iteration at t0. */
    /* rhs = hN*f0 */
    for (i = 0; i < nx; i++) {
      Delta[i] = hN*f0[i];
    }

    /* Delta = (U \ (L \ rhs)) */
    rt_ForwardSubstitutionRR_Dbl(W, Delta,
      f1, nx,
      1, pivots,
      1);
    rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
      Delta, nx,
      1, 0);

    /* ytmp = y0 + Delta
       ztmp = (ytmp-y0)/h
     */
    (void) memcpy(x1, x0,
                  (uint_T)nx*sizeof(real_T));
    for (i = 0; i < nx; i++) {
      x1[i] += Delta[i];
      ztmp[i] = Delta[i]/hN;
    }

    /* Additional Newton's iterations, if desired.
       for iter = 2:NewtIter
       rhs = hN*feval(odefun,tn,ytmp,extraArgs{:}) - M*(ytmp - yn);
       if statedepM   % only for state dep. Mdel ~= 0
       Mdel = M - feval(massfun,tn,ytmp);
       rhs = rhs + Mdel*ztmp*h;
       end
       Delta = ( U \ ( L \ rhs ) );
       ytmp = ytmp + Delta;
       ztmp = (ytmp - yn)/h
       end
     */
    rtsiSetT(si, t0);
    rtsiSetdX(si, f1);
    for (iter = 1; iter < numIter; iter++) {
      Subsystem_step();
      Subsystem_derivatives();
      for (i = 0; i < nx; i++) {
        Delta[i] = hN*f1[i];
        xtmp[i] = x1[i] - x0[i];
      }

      /* rhs = hN*f(tn,ytmp) - (ytmp-yn) */
      for (i = 0; i < nx; i++) {
        Delta[i] -= xtmp[i];
      }

      rt_ForwardSubstitutionRR_Dbl(W, Delta,
        f1, nx,
        1, pivots,
        1);
      rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
        Delta, nx,
        1, 0);

      /* ytmp = ytmp + delta
         ztmp = (ytmp - yn)/h
       */
      for (i = 0; i < nx; i++) {
        x1[i] += Delta[i];
        ztmp[i] = (x1[i] - x0[i])/hN;
      }
    }

    /* Steps from t0+hN to t1 -- subintegration of N(j) steps for extrapolation
       ttmp = t0;
       for i = 2:N(j)
       ttmp = ttmp + hN
       ytmp0 = ytmp;
       for iter = 1:NewtIter
       rhs = (ytmp0 - ytmp) + hN*feval(odefun,ttmp,ytmp,extraArgs{:});
       Delta = ( U \ ( L \ rhs ) );
       ytmp = ytmp + Delta;
       end
       end
     */
    for (k = 1; k < N[j]; k++) {
      t1 = t0 + k*hN;
      (void) memcpy(x1start, x1,
                    (uint_T)nx*sizeof(real_T));
      rtsiSetT(si, t1);
      rtsiSetdX(si, f1);
      for (iter = 0; iter < numIter; iter++) {
        Subsystem_step();
        Subsystem_derivatives();
        if (iter == 0) {
          for (i = 0; i < nx; i++) {
            Delta[i] = hN*f1[i];
          }
        } else {
          for (i = 0; i < nx; i++) {
            Delta[i] = hN*f1[i];
            xtmp[i] = (x1[i]-x1start[i]);
          }

          /* rhs = hN*f(tn,ytmp) - M*(ytmp-yn) */
          for (i = 0; i < nx; i++) {
            Delta[i] -= xtmp[i];
          }
        }

        rt_ForwardSubstitutionRR_Dbl(W, Delta,
          f1, nx,
          1, pivots,
          1);
        rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
          Delta, nx,
          1, 0);

        /* ytmp = ytmp + Delta
           ztmp = (ytmp - ytmp0)/h
         */
        for (i = 0; i < nx; i++) {
          x1[i] += Delta[i];
          ztmp[i] = (x1[i] - x1start[i])/hN;
        }
      }
    }

    /* Extrapolate to order j
       E(:,j) = ytmp
       for k = j:-1:2
       coef = N(k-1)/(N(j) - N(k-1))
       E(:,k-1) = E(:,k) + coef*( E(:,k) - E(:,k-1) )
       end
     */
    (void) memcpy(&(E[nx*j]), x1,
                  (uint_T)nx*sizeof(real_T));
    for (k = j; k > 0; k--) {
      real_T coef = (real_T)(N[k-1]) / (N[j]-N[k-1]);
      for (i = 0; i < nx; i++) {
        x1[i] = E[nx*k+i] + coef*(E[nx*k+i] - E[nx*(k-1)+i]);
      }

      (void) memcpy(&(E[nx*(k-1)]), x1,
                    (uint_T)nx*sizeof(real_T));
    }
  }

  /* x1 = E(:,1); */
  (void) memcpy(x1, E,
                (uint_T)nx*sizeof(real_T));

  /* t1 = t0 + h; */
  rtsiSetT(si,rtsiGetSolverStopTime(si));
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    real_T tmp;
    real_T tmp_0;
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* Model step function */
void Subsystem_step(void)
{
  real_T Divide1_j;
  real_T Product_o;
  real_T Switch_idx_0;
  int32_T i1;
  int32_T k1;
  int32_T k2;
  int32_T kq;
  static const real_T b[102] = { 0.0, 381.82, 534.53, 686.8, 763.6, 836.5, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.2, 1.2, 1.2, 1.2, 1.2, 1.0, 1.55, 1.5328,
    1.5321, 1.5316, 1.5327, 1.0, 1.7, 1.73, 1.7418, 1.7413, 1.7409, 1.0, 2.0216,
    2.0172, 2.017, 2.0169, 2.0168, 1.0, 3.4511, 3.4473, 3.4496, 3.4509, 3.4509,
    1.0, 4.7626, 4.753, 4.7589, 4.7605, 4.8437, 1.0, 6.1721, 6.2159, 6.3032,
    6.3527, 6.344, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.22, 0.225, 0.22, 0.215,
    0.216, 0.0, 0.292801, 0.296743, 0.2975, 0.2977, 0.2987, 0.0, 0.30376,
    0.306989, 0.3091, 0.309585, 0.310585, 0.0, 0.311119, 0.315387, 0.316097,
    0.3163, 0.3173, 0.0, 0.318322, 0.320272, 0.320966, 0.321209, 0.322209, 0.0,
    0.31856, 0.320411, 0.321058, 0.321212, 0.322212, 0.0, 0.318617, 0.320429,
    0.321185, 0.321266, 0.322266 };

  static const real_T b_0[120] = { 0.0, 568.67, 853.0, 995.78, 1138.0, 1279.5,
    1422.6, 1558.3, 1.0, 1.0, 1.22, 1.36, 1.54, 1.8, 2.0, 2.15, 1.0, 1.0, 1.31,
    1.49, 1.79, 2.5, 3.02, 4.25, 1.0, 1.09, 1.4, 1.69, 1.99, 2.83, 3.37, 4.71,
    1.0, 1.09, 1.5, 1.89, 2.3, 3.04, 3.82, 4.82, 1.0, 1.18, 1.6, 2.0, 2.46, 3.17,
    3.96, 4.93, 1.0, 1.23, 1.7, 2.08, 2.53, 3.22, 4.02, 4.96, 1.0, 1.28, 1.78,
    2.12, 2.59, 3.22, 4.02, 4.96, 0.0, 0.32, 0.4328, 0.5063, 0.5899, 0.677,
    0.7482, 0.7916, 0.0, 0.32, 0.4285, 0.5063, 0.5898, 0.6702, 0.7448, 0.7843,
    0.0, 0.2918, 0.4194, 0.4956, 0.5848, 0.6464, 0.736, 0.7542, 0.0, 0.2918,
    0.4049, 0.4649, 0.5524, 0.5914, 0.6792, 0.7316, 0.0, 0.2613, 0.3803, 0.429,
    0.507, 0.5265, 0.6247, 0.6913, 0.0, 0.2295, 0.3376, 0.3899, 0.4679, 0.4675,
    0.5848, 0.6616, 0.0, 0.1977, 0.264, 0.3112, 0.3478, 0.3897, 0.4294, 0.5072 };

  static const real_T b_1[85] = { 381.82, 534.53, 686.8, 763.6, 836.5, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.12, 1.2, 1.3, 1.5316, 1.5327, 1.5334, 1.5328, 1.6, 1.7742,
    1.8663, 1.7, 1.7419, 2.017, 2.0169, 2.2, 2.0216, 2.0172, 2.7374, 2.7373,
    2.8255, 3.4511, 3.4473, 3.4496, 3.4509, 3.4509, 4.7626, 4.753, 4.7589,
    4.7605, 4.8437, 6.1721, 6.2159, 6.3032, 6.3527, 6.344, 0.5, 0.5, 0.5, 0.5,
    0.5, 0.72, 0.74, 0.74, 0.7387, 0.7, 0.6487, 0.744, 0.78, 0.7591, 0.74, 0.61,
    0.7212, 0.7659, 0.7794, 0.78, 0.5727, 0.6947, 0.7333, 0.7596, 0.7574, 0.4734,
    0.5992, 0.6849, 0.7129, 0.7349, 0.4179, 0.5305, 0.6093, 0.6376, 0.603,
    0.3758, 0.476, 0.5453, 0.5599, 0.3884 };

  static const real_T b_2[77] = { 568.67, 853.0, 995.78, 1138.0, 1279.5, 1422.6,
    1558.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1977, 0.3376, 0.3112, 0.4679,
    0.5079, 0.6247, 0.6913, 0.2613, 0.4048, 0.429, 0.5524, 0.6464, 0.7283,
    0.7787, 0.2918, 0.4196, 0.5063, 0.5899, 0.6758, 0.7469, 0.7915, 0.32, 0.4328,
    0.507, 0.59, 0.677, 0.749, 0.792, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75,
    0.8264, 0.77, 0.7946, 0.7655, 0.7515, 0.7241, 0.8289, 0.7396, 0.8218, 0.7963,
    0.7846, 0.7491, 0.7191, 0.7069, 0.6439, 0.561, 0.6148, 0.6969, 0.6977,
    0.6228, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4 };

  boolean_T guard1 = false;
  if (rtmIsMajorTimeStep(Subsystem_M)) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&Subsystem_M->solverInfo,
                          ((Subsystem_M->Timing.clockTick0+1)*
      Subsystem_M->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(Subsystem_M)) {
    Subsystem_M->Timing.t[0] = rtsiGetT(&Subsystem_M->solverInfo);
  }

  /* MinMax: '<S4>/Min1' incorporates:
   *  Constant: '<S4>/Constant'
   *  Constant: '<S4>/n_max'
   *  Inport: '<Root>/G_fuel'
   *  MinMax: '<S4>/Min'
   */
  Subsystem_B.Min1 = fmax(fmin(Subsystem_P.n_max_Value, Subsystem_U.G_fuel),
    Subsystem_P.Constant_Value);
  if (rtmIsMajorTimeStep(Subsystem_M)) {
    /* Sum: '<S2>/Sum' incorporates:
     *  Constant: '<S1>/������ (�)'
     *  Constant: '<S2>/Constant'
     *  Gain: '<S2>/Gain1'
     */
    Subsystem_B.T = Subsystem_P.Constant_Value_g - Subsystem_P.Gain1_Gain *
      Subsystem_P._Value;

    /* Sum: '<S2>/Sum1' incorporates:
     *  Constant: '<S1>/������ (�)'
     *  Constant: '<S2>/Constant1'
     *  Gain: '<S2>/Gain2'
     */
    Subsystem_B.MathFunction = Subsystem_P.Constant1_Value -
      Subsystem_P.Gain2_Gain * Subsystem_P._Value;

    /* Math: '<S2>/Math Function' incorporates:
     *  Constant: '<S2>/Constant2'
     */
    if ((Subsystem_B.MathFunction < 0.0) && (Subsystem_P.Constant2_Value > floor
         (Subsystem_P.Constant2_Value))) {
      Subsystem_B.MathFunction = -rt_powd_snf(-Subsystem_B.MathFunction,
        Subsystem_P.Constant2_Value);
    } else {
      Subsystem_B.MathFunction = rt_powd_snf(Subsystem_B.MathFunction,
        Subsystem_P.Constant2_Value);
    }

    /* End of Math: '<S2>/Math Function' */

    /* Gain: '<S2>/Gain' */
    Subsystem_B.p = Subsystem_P.Gain_Gain * Subsystem_B.MathFunction;
  }

  /* Integrator: '<S12>/Integrator' */
  if (Subsystem_DW.Integrator_IWORK != 0) {
    Subsystem_X.Integrator_CSTATE = Subsystem_B.T;
  }

  /* Integrator: '<S11>/Integrator' */
  if (Subsystem_DW.Integrator_IWORK_a != 0) {
    Subsystem_X.Integrator_CSTATE_e = Subsystem_B.p;
  }

  /* Fcn: '<S10>/Fcn1' incorporates:
   *  Integrator: '<S12>/Integrator'
   */
  Subsystem_B.MathFunction = 288.0 / Subsystem_X.Integrator_CSTATE;
  if (Subsystem_B.MathFunction < 0.0) {
    Subsystem_B.MathFunction = -rt_powd_snf(-Subsystem_B.MathFunction, 0.5);
  } else {
    Subsystem_B.MathFunction = rt_powd_snf(Subsystem_B.MathFunction, 0.5);
  }

  if (rtmIsMajorTimeStep(Subsystem_M)) {
    /* Gain: '<S9>/Gain' */
    Subsystem_B.Gain = 1.0 / Subsystem_P.sig_s * Subsystem_B.p;
  }

  /* Product: '<S10>/Divide' incorporates:
   *  Integrator: '<S11>/Integrator'
   */
  Subsystem_B.Divide = 1.0 / Subsystem_B.Gain * Subsystem_X.Integrator_CSTATE_e;

  /* ManualSwitch: '<S8>/Manual Switch' incorporates:
   *  Integrator: '<S8>/Integrator'
   *  Saturate: '<S8>/Saturation1'
   */
  if (Subsystem_P.ManualSwitch_CurrentSetting == 1) {
    /* Saturate: '<S8>/Saturation' incorporates:
     *  Integrator: '<S8>/Integrator'
     */
    if (Subsystem_X.Integrator_CSTATE_o > Subsystem_P.Saturation_UpperSat) {
      /* ManualSwitch: '<S8>/Manual Switch' */
      Subsystem_B.ManualSwitch = Subsystem_P.Saturation_UpperSat;
    } else if (Subsystem_X.Integrator_CSTATE_o < Subsystem_P.Saturation_LowerSat)
    {
      /* ManualSwitch: '<S8>/Manual Switch' */
      Subsystem_B.ManualSwitch = Subsystem_P.Saturation_LowerSat;
    } else {
      /* ManualSwitch: '<S8>/Manual Switch' */
      Subsystem_B.ManualSwitch = Subsystem_X.Integrator_CSTATE_o;
    }

    /* End of Saturate: '<S8>/Saturation' */
  } else if (Subsystem_X.Integrator_CSTATE_o > Subsystem_P.Saturation1_UpperSat)
  {
    /* Saturate: '<S8>/Saturation1' incorporates:
     *  ManualSwitch: '<S8>/Manual Switch'
     */
    Subsystem_B.ManualSwitch = Subsystem_P.Saturation1_UpperSat;
  } else if (Subsystem_X.Integrator_CSTATE_o < Subsystem_P.Saturation1_LowerSat)
  {
    /* Saturate: '<S8>/Saturation1' incorporates:
     *  ManualSwitch: '<S8>/Manual Switch'
     */
    Subsystem_B.ManualSwitch = Subsystem_P.Saturation1_LowerSat;
  } else {
    /* ManualSwitch: '<S8>/Manual Switch' incorporates:
     *  Integrator: '<S8>/Integrator'
     *  Saturate: '<S8>/Saturation1'
     */
    Subsystem_B.ManualSwitch = Subsystem_X.Integrator_CSTATE_o;
  }

  /* End of ManualSwitch: '<S8>/Manual Switch' */

  /* Gain: '<S8>/Gain' */
  Subsystem_Y.n = Subsystem_P.Gain_Gain_h * Subsystem_B.ManualSwitch;

  /* Product: '<S10>/Product' incorporates:
   *  Fcn: '<S10>/Fcn1'
   */
  Subsystem_B.Product = Subsystem_B.MathFunction * Subsystem_Y.n;

  /* MATLAB Function: '<S10>/��������������  �������' */
  k1 = -1;
  k2 = -1;
  Subsystem_B.XA = 0.0;
  Subsystem_B.XB = 0.0;
  Subsystem_B.G_pr = 0.0;
  if (Subsystem_B.Product > 836.5) {
    Subsystem_B.XA = 763.6;
    Subsystem_B.XB = 836.5;
    k1 = 10;
    k2 = 11;
  }

  if (Subsystem_B.Product < 0.0) {
    Subsystem_B.XA = 0.0;
    Subsystem_B.XB = 381.82;
    k1 = 6;
    k2 = 7;
    for (i1 = 0; i1 < 16; i1++) {
      kq = i1 * 6 + 6;
      Subsystem_B.Gkpd[i1] = (b[i1 * 6 + 7] - b[kq]) * Subsystem_B.Product /
        381.82 + b[kq];
    }
  } else {
    for (kq = 0; kq < 5; kq++) {
      if (Subsystem_B.Product >= b[kq]) {
        Switch_idx_0 = b[kq + 1];
        if (Subsystem_B.Product <= Switch_idx_0) {
          Subsystem_B.XA = b[kq];
          Subsystem_B.XB = Switch_idx_0;
          k1 = kq + 6;
          k2 = kq + 7;
        }
      }
    }

    for (i1 = 0; i1 < 16; i1++) {
      kq = i1 * 6 + k1;
      if (Subsystem_B.XB == Subsystem_B.XA) {
        Subsystem_B.XB += 0.1;
      }

      Subsystem_B.Gkpd[i1] = (b[i1 * 6 + k2] - b[kq]) * (Subsystem_B.Product -
        Subsystem_B.XA) / (Subsystem_B.XB - Subsystem_B.XA) + b[kq];
    }
  }

  if (Subsystem_B.Divide > Subsystem_B.Gkpd[7]) {
    Subsystem_B.XA = Subsystem_B.Gkpd[6];
    Subsystem_B.XB = Subsystem_B.Gkpd[7];
    k1 = 14;
    k2 = 15;
    if (Subsystem_B.Gkpd[6] == Subsystem_B.Gkpd[7]) {
      Subsystem_B.XB = Subsystem_B.Gkpd[7] + 0.1;
    }

    Subsystem_B.G_pr = (Subsystem_B.Gkpd[15] - Subsystem_B.Gkpd[14]) *
      (Subsystem_B.Divide - Subsystem_B.Gkpd[6]) / (Subsystem_B.XB -
      Subsystem_B.Gkpd[6]) + Subsystem_B.Gkpd[14];
  }

  if (Subsystem_B.Divide < Subsystem_B.Gkpd[0]) {
    Subsystem_B.XB = Subsystem_B.Gkpd[1];
    if (Subsystem_B.Gkpd[0] == Subsystem_B.Gkpd[1]) {
      Subsystem_B.XB = Subsystem_B.Gkpd[1] + 0.1;
    }

    Subsystem_B.G_pr = (Subsystem_B.Gkpd[9] - Subsystem_B.Gkpd[8]) *
      (Subsystem_B.Divide - Subsystem_B.Gkpd[0]) / (Subsystem_B.XB -
      Subsystem_B.Gkpd[0]) + Subsystem_B.Gkpd[8];
  } else {
    for (kq = 0; kq < 7; kq++) {
      guard1 = false;
      if (Subsystem_B.Divide >= Subsystem_B.Gkpd[kq]) {
        Switch_idx_0 = Subsystem_B.Gkpd[kq + 1];
        if (Subsystem_B.Divide <= Switch_idx_0) {
          Subsystem_B.XA = Subsystem_B.Gkpd[kq];
          Subsystem_B.XB = Switch_idx_0;
          k1 = kq + 8;
          k2 = kq + 9;
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        if (Subsystem_B.XB == Subsystem_B.XA) {
          Subsystem_B.XB += 0.1;
        }

        Subsystem_B.G_pr = (Subsystem_B.Gkpd[k2] - Subsystem_B.Gkpd[k1]) *
          (Subsystem_B.Divide - Subsystem_B.XA) / (Subsystem_B.XB -
          Subsystem_B.XA) + Subsystem_B.Gkpd[k1];
      }
    }
  }

  /* End of MATLAB Function: '<S10>/��������������  �������' */

  /* Product: '<S10>/Product1' incorporates:
   *  Fcn: '<S10>/Fcn1'
   *  Gain: '<S10>/Gain2'
   *  Integrator: '<S11>/Integrator'
   */
  Subsystem_B.G_pr = Subsystem_P.Gain2_Gain_j * Subsystem_X.Integrator_CSTATE_e *
    Subsystem_B.G_pr * Subsystem_B.MathFunction;

  /* Switch: '<S15>/Switch' */
  if (Subsystem_B.Min1 > Subsystem_P.Switch_Threshold) {
    /* Product: '<S15>/Divide1' incorporates:
     *  Integrator: '<S12>/Integrator'
     *  Lookup_n-D: '<S15>/i_g'
     */
    Divide1_j = look1_binlxpw(Subsystem_X.Integrator_CSTATE, Subsystem_P.T_g_arr,
      Subsystem_P.i_g_arr, 7U) / Subsystem_X.Integrator_CSTATE;

    /* Switch: '<S15>/Switch' incorporates:
     *  Integrator: '<S12>/Integrator'
     *  Product: '<S15>/Product'
     */
    Switch_idx_0 = Subsystem_B.G_pr * Divide1_j * Subsystem_X.Integrator_CSTATE;
  } else {
    /* Product: '<S15>/Divide2' incorporates:
     *  Integrator: '<S12>/Integrator'
     *  Lookup_n-D: '<S15>/i_k'
     */
    Divide1_j = look1_binlxpw(Subsystem_X.Integrator_CSTATE, Subsystem_P.T_k_arr,
      Subsystem_P.i_k_arr, 7U) / Subsystem_X.Integrator_CSTATE;

    /* Switch: '<S15>/Switch' incorporates:
     *  Integrator: '<S12>/Integrator'
     *  Product: '<S15>/Product1'
     */
    Switch_idx_0 = Subsystem_B.G_pr * Divide1_j * Subsystem_X.Integrator_CSTATE;
  }

  /* End of Switch: '<S15>/Switch' */
  if (rtmIsMajorTimeStep(Subsystem_M)) {
    /* Fcn: '<S5>/Fcn' incorporates:
     *  Constant: '<S1>/�������� (���)'
     *  Fcn: '<S5>/Fcn1'
     */
    Subsystem_B.G_k = 0.19999999999999996 * rt_powd_snf(Subsystem_P._Value_a,
      2.0) + 1.0;

    /* Product: '<S5>/Product' incorporates:
     *  Fcn: '<S5>/Fcn'
     */
    Subsystem_B.T_vh = Subsystem_B.G_k * Subsystem_B.T;

    /* Fcn: '<S7>/Fcn' */
    Subsystem_B.MathFunction = 288.0 / Subsystem_B.T_vh;
    if (Subsystem_B.MathFunction < 0.0) {
      /* Fcn: '<S7>/Fcn' */
      Subsystem_B.Fcn = -sqrt(-Subsystem_B.MathFunction);
    } else {
      /* Fcn: '<S7>/Fcn' */
      Subsystem_B.Fcn = sqrt(Subsystem_B.MathFunction);
    }

    /* End of Fcn: '<S7>/Fcn' */

    /* Fcn: '<S5>/Fcn1' */
    if (Subsystem_B.G_k < 0.0) {
      Subsystem_B.MathFunction = -rt_powd_snf(-Subsystem_B.G_k,
        3.5000000000000004);
    } else {
      Subsystem_B.MathFunction = rt_powd_snf(Subsystem_B.G_k, 3.5000000000000004);
    }

    /* Gain: '<S5>/Gain1' incorporates:
     *  Product: '<S5>/Product1'
     */
    Subsystem_B.p_vh = Subsystem_B.MathFunction * Subsystem_B.p *
      Subsystem_P.sig_vh;
  }

  /* Product: '<S7>/Divide' incorporates:
   *  Gain: '<S11>/Gain1'
   *  Integrator: '<S11>/Integrator'
   */
  Subsystem_B.Divide_o = 1.0 / Subsystem_P.sig_g *
    Subsystem_X.Integrator_CSTATE_e / Subsystem_B.p_vh;

  /* Product: '<S7>/Product' */
  Product_o = Subsystem_Y.n * Subsystem_B.Fcn;

  /* MATLAB Function: '<S7>/��������������  �����������' */
  k1 = -1;
  k2 = -1;
  Subsystem_B.XA = 0.0;
  Subsystem_B.XB = 0.0;
  Subsystem_B.G_pr_l = 0.0;
  if (Product_o > 1558.3) {
    Subsystem_B.XA = 1422.6;
    Subsystem_B.XB = 1558.3;
    k1 = 14;
    k2 = 15;
  }

  if (Product_o < 0.0) {
    Subsystem_B.XA = 0.0;
    Subsystem_B.XB = 568.67;
    k1 = 8;
    k2 = 9;
    for (i1 = 0; i1 < 14; i1++) {
      int32_T kq_tmp;
      kq_tmp = i1 << 3;
      Subsystem_B.G_k = b_0[kq_tmp + 8];
      Subsystem_B.Gkpd_m[i1] = (b_0[kq_tmp + 9] - Subsystem_B.G_k) * Product_o /
        568.67 + Subsystem_B.G_k;
    }
  } else {
    for (kq = 0; kq < 7; kq++) {
      if (Product_o >= b_0[kq]) {
        Subsystem_B.MathFunction = b_0[kq + 1];
        if (Product_o <= Subsystem_B.MathFunction) {
          Subsystem_B.XA = b_0[kq];
          Subsystem_B.XB = Subsystem_B.MathFunction;
          k1 = kq + 8;
          k2 = kq + 9;
        }
      }
    }

    for (i1 = 0; i1 < 14; i1++) {
      int32_T kq_tmp;
      kq_tmp = i1 << 3;
      kq = kq_tmp + k1;
      if (Subsystem_B.XB == Subsystem_B.XA) {
        Subsystem_B.XB += 0.1;
      }

      Subsystem_B.Gkpd_m[i1] = (b_0[kq_tmp + k2] - b_0[kq]) * (Product_o -
        Subsystem_B.XA) / (Subsystem_B.XB - Subsystem_B.XA) + b_0[kq];
    }
  }

  if (Subsystem_B.Divide_o > Subsystem_B.Gkpd_m[6]) {
    Subsystem_B.XA = Subsystem_B.Gkpd_m[5];
    Subsystem_B.XB = Subsystem_B.Gkpd_m[6];
    k1 = 12;
    k2 = 13;
    if (Subsystem_B.Gkpd_m[5] == Subsystem_B.Gkpd_m[6]) {
      Subsystem_B.XB = Subsystem_B.Gkpd_m[6] + 0.1;
    }

    Subsystem_B.G_pr_l = (Subsystem_B.Gkpd_m[13] - Subsystem_B.Gkpd_m[12]) *
      (Subsystem_B.Divide_o - Subsystem_B.Gkpd_m[5]) / (Subsystem_B.XB -
      Subsystem_B.Gkpd_m[5]) + Subsystem_B.Gkpd_m[12];
  }

  if (Subsystem_B.Divide_o < Subsystem_B.Gkpd_m[0]) {
    Subsystem_B.XB = Subsystem_B.Gkpd_m[1];
    if (Subsystem_B.Gkpd_m[0] == Subsystem_B.Gkpd_m[1]) {
      Subsystem_B.XB = Subsystem_B.Gkpd_m[1] + 0.1;
    }

    Subsystem_B.G_pr_l = (Subsystem_B.Gkpd_m[8] - Subsystem_B.Gkpd_m[7]) *
      (Subsystem_B.Divide_o - Subsystem_B.Gkpd_m[0]) / (Subsystem_B.XB -
      Subsystem_B.Gkpd_m[0]) + Subsystem_B.Gkpd_m[7];
  } else {
    for (kq = 0; kq < 6; kq++) {
      guard1 = false;
      if (Subsystem_B.Divide_o >= Subsystem_B.Gkpd_m[kq]) {
        Subsystem_B.MathFunction = Subsystem_B.Gkpd_m[kq + 1];
        if (Subsystem_B.Divide_o <= Subsystem_B.MathFunction) {
          Subsystem_B.XA = Subsystem_B.Gkpd_m[kq];
          Subsystem_B.XB = Subsystem_B.MathFunction;
          k1 = kq + 7;
          k2 = kq + 8;
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        if (Subsystem_B.XB == Subsystem_B.XA) {
          Subsystem_B.XB += 0.1;
        }

        Subsystem_B.G_pr_l = (Subsystem_B.Gkpd_m[k2] - Subsystem_B.Gkpd_m[k1]) *
          (Subsystem_B.Divide_o - Subsystem_B.XA) / (Subsystem_B.XB -
          Subsystem_B.XA) + Subsystem_B.Gkpd_m[k1];
      }
    }
  }

  /* End of MATLAB Function: '<S7>/��������������  �����������' */

  /* Saturate: '<S7>/Saturation' */
  if (Subsystem_B.G_pr_l > Subsystem_P.Saturation_UpperSat_k) {
    /* Saturate: '<S7>/Saturation' */
    Subsystem_B.G_pr_l = Subsystem_P.Saturation_UpperSat_k;
  } else if (Subsystem_B.G_pr_l < Subsystem_P.Saturation_LowerSat_b) {
    /* Saturate: '<S7>/Saturation' */
    Subsystem_B.G_pr_l = Subsystem_P.Saturation_LowerSat_b;
  }

  /* End of Saturate: '<S7>/Saturation' */
  if (rtmIsMajorTimeStep(Subsystem_M)) {
    /* Gain: '<S7>/Gain1' */
    Subsystem_B.Gain1_m = Subsystem_P.Gain1_Gain_b * Subsystem_B.p_vh;
  }

  /* Product: '<S7>/Product1' */
  Subsystem_B.G_k = Subsystem_B.G_pr_l * Subsystem_B.Fcn * Subsystem_B.Gain1_m;

  /* Lookup_n-D: '<S6>/R' incorporates:
   *  Product: '<S6>/K'
   */
  Subsystem_B.R = look1_binlxpw(Subsystem_B.G_k / Subsystem_B.Min1,
    Subsystem_P.K_arr, Subsystem_P.R_arr, 11U);

  /* Product: '<S6>/Divide' incorporates:
   *  Sum: '<S6>/Sum'
   */
  Subsystem_B.Divide_b = Divide1_j / (Divide1_j - Subsystem_B.R);

  /* Fcn: '<S10>/Fcn3' */
  Subsystem_B.uDLookupTable = (Subsystem_B.Divide_b - 1.0) /
    Subsystem_B.Divide_b;

  /* Math: '<S10>/Math Function1'
   *
   * About '<S10>/Math Function1':
   *  Operator: reciprocal
   */
  Subsystem_B.i_k_p = 1.0 / Subsystem_B.Divide;

  /* MATLAB Function: '<S10>/��� �������' */
  k1 = -1;
  k2 = -1;
  Subsystem_B.XA = 0.0;
  Subsystem_B.XB = 0.0;
  Subsystem_B.kpd = 0.0;
  if (Subsystem_B.Product > 836.5) {
    Subsystem_B.XA = 763.6;
    Subsystem_B.XB = 836.5;
    k1 = 8;
    k2 = 9;
  }

  if (Subsystem_B.Product < 381.82) {
    Subsystem_B.XA = 381.82;
    Subsystem_B.XB = 534.53;
    k1 = 5;
    k2 = 6;
    for (i1 = 0; i1 < 16; i1++) {
      kq = i1 * 5 + 5;
      Subsystem_B.Gkpd[i1] = (b_1[i1 * 5 + 6] - b_1[kq]) * (Subsystem_B.Product
        - 381.82) / 152.70999999999998 + b_1[kq];
    }
  } else {
    if ((Subsystem_B.Product >= 381.82) && (Subsystem_B.Product <= 534.53)) {
      Subsystem_B.XA = 381.82;
      Subsystem_B.XB = 534.53;
      k1 = 5;
      k2 = 6;
    }

    if ((Subsystem_B.Product >= 534.53) && (Subsystem_B.Product <= 686.8)) {
      Subsystem_B.XA = 534.53;
      Subsystem_B.XB = 686.8;
      k1 = 6;
      k2 = 7;
    }

    if ((Subsystem_B.Product >= 686.8) && (Subsystem_B.Product <= 763.6)) {
      Subsystem_B.XA = 686.8;
      Subsystem_B.XB = 763.6;
      k1 = 7;
      k2 = 8;
    }

    if ((Subsystem_B.Product >= 763.6) && (Subsystem_B.Product <= 836.5)) {
      Subsystem_B.XA = 763.6;
      Subsystem_B.XB = 836.5;
      k1 = 8;
      k2 = 9;
    }

    for (i1 = 0; i1 < 16; i1++) {
      kq = i1 * 5 + k1;
      if (Subsystem_B.XB == Subsystem_B.XA) {
        Subsystem_B.XB += 0.1;
      }

      Subsystem_B.Gkpd[i1] = (b_1[i1 * 5 + k2] - b_1[kq]) * (Subsystem_B.Product
        - Subsystem_B.XA) / (Subsystem_B.XB - Subsystem_B.XA) + b_1[kq];
    }
  }

  if (Subsystem_B.Divide > Subsystem_B.Gkpd[7]) {
    Subsystem_B.XA = Subsystem_B.Gkpd[6];
    Subsystem_B.XB = Subsystem_B.Gkpd[7];
    k1 = 14;
    k2 = 15;
    if (Subsystem_B.Gkpd[6] == Subsystem_B.Gkpd[7]) {
      Subsystem_B.XB = Subsystem_B.Gkpd[7] + 0.1;
    }

    Subsystem_B.kpd = (Subsystem_B.Gkpd[15] - Subsystem_B.Gkpd[14]) *
      (Subsystem_B.Divide - Subsystem_B.Gkpd[6]) / (Subsystem_B.XB -
      Subsystem_B.Gkpd[6]) + Subsystem_B.Gkpd[14];
  }

  if (Subsystem_B.Divide < Subsystem_B.Gkpd[0]) {
    Subsystem_B.XB = Subsystem_B.Gkpd[1];
    if (Subsystem_B.Gkpd[0] == Subsystem_B.Gkpd[1]) {
      Subsystem_B.XB = Subsystem_B.Gkpd[1] + 0.1;
    }

    Subsystem_B.kpd = (Subsystem_B.Gkpd[9] - Subsystem_B.Gkpd[8]) *
      (Subsystem_B.Divide - Subsystem_B.Gkpd[0]) / (Subsystem_B.XB -
      Subsystem_B.Gkpd[0]) + Subsystem_B.Gkpd[8];
  } else {
    for (kq = 0; kq < 7; kq++) {
      guard1 = false;
      if (Subsystem_B.Divide >= Subsystem_B.Gkpd[kq]) {
        Subsystem_B.MathFunction = Subsystem_B.Gkpd[kq + 1];
        if (Subsystem_B.Divide <= Subsystem_B.MathFunction) {
          Subsystem_B.XA = Subsystem_B.Gkpd[kq];
          Subsystem_B.XB = Subsystem_B.MathFunction;
          k1 = kq + 8;
          k2 = kq + 9;
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        if (Subsystem_B.XB == Subsystem_B.XA) {
          Subsystem_B.XB += 0.1;
        }

        Subsystem_B.kpd = (Subsystem_B.Gkpd[k2] - Subsystem_B.Gkpd[k1]) *
          (Subsystem_B.Divide - Subsystem_B.XA) / (Subsystem_B.XB -
          Subsystem_B.XA) + Subsystem_B.Gkpd[k1];
      }
    }
  }

  /* End of MATLAB Function: '<S10>/��� �������' */

  /* Math: '<S10>/Math Function2' */
  if ((Subsystem_B.i_k_p < 0.0) && (Subsystem_B.uDLookupTable > floor
       (Subsystem_B.uDLookupTable))) {
    Subsystem_B.i_k_p = -rt_powd_snf(-Subsystem_B.i_k_p,
      Subsystem_B.uDLookupTable);
  } else {
    Subsystem_B.i_k_p = rt_powd_snf(Subsystem_B.i_k_p, Subsystem_B.uDLookupTable);
  }

  /* End of Math: '<S10>/Math Function2' */

  /* Product: '<S10>/Product4' incorporates:
   *  Fcn: '<S10>/Fcn5'
   *  Fcn: '<S10>/Fcn6'
   *  Integrator: '<S12>/Integrator'
   *  Product: '<S10>/Product3'
   */
  Subsystem_B.Divide = (1.0 - (1.0 - Subsystem_B.i_k_p) * Subsystem_B.kpd) *
    Subsystem_X.Integrator_CSTATE;

  /* Outport: '<Root>/P' incorporates:
   *  Gain: '<S9>/Gain1'
   *  Gain: '<S9>/Gain2'
   *  Product: '<S9>/Divide'
   *  Product: '<S9>/Divide1'
   *  Product: '<S9>/Product'
   */
  Subsystem_Y.P = Subsystem_B.G_pr / (Subsystem_B.Gain / Subsystem_B.R /
    Subsystem_B.Divide * Subsystem_P.F_s) * Subsystem_P.fi * Subsystem_B.G_pr;

  /* MATLAB Function: '<S7>/��� �����������' */
  k1 = -1;
  k2 = -1;
  Subsystem_B.XA = 0.0;
  Subsystem_B.XB = 0.0;
  Subsystem_B.Product = 0.0;
  if (Product_o > 1558.3) {
    Subsystem_B.XA = 1422.6;
    Subsystem_B.XB = 1558.3;
    k1 = 12;
    k2 = 13;
  }

  if (Product_o < 568.67) {
    Subsystem_B.XA = 568.67;
    Subsystem_B.XB = 853.0;
    k1 = 7;
    k2 = 8;
    for (i1 = 0; i1 < 10; i1++) {
      kq = i1 * 7 + 7;
      Subsystem_B.Gkpd_c[i1] = (b_2[i1 * 7 + 8] - b_2[kq]) * (Product_o - 568.67)
        / 284.33000000000004 + b_2[kq];
    }
  } else {
    for (kq = 0; kq < 6; kq++) {
      if (Product_o >= b_2[kq]) {
        Subsystem_B.MathFunction = b_2[kq + 1];
        if (Product_o <= Subsystem_B.MathFunction) {
          Subsystem_B.XA = b_2[kq];
          Subsystem_B.XB = Subsystem_B.MathFunction;
          k1 = kq + 7;
          k2 = kq + 8;
        }
      }
    }

    for (i1 = 0; i1 < 10; i1++) {
      kq = i1 * 7 + k1;
      if (Subsystem_B.XB == Subsystem_B.XA) {
        Subsystem_B.XB += 0.1;
      }

      Subsystem_B.Gkpd_c[i1] = (b_2[i1 * 7 + k2] - b_2[kq]) * (Product_o -
        Subsystem_B.XA) / (Subsystem_B.XB - Subsystem_B.XA) + b_2[kq];
    }
  }

  if (Subsystem_B.G_pr_l > Subsystem_B.Gkpd_c[4]) {
    Subsystem_B.XA = Subsystem_B.Gkpd_c[3];
    Subsystem_B.XB = Subsystem_B.Gkpd_c[4];
    k1 = 8;
    k2 = 9;
    if (Subsystem_B.Gkpd_c[3] == Subsystem_B.Gkpd_c[4]) {
      Subsystem_B.XB = Subsystem_B.Gkpd_c[4] + 0.1;
    }

    Subsystem_B.Product = (Subsystem_B.Gkpd_c[9] - Subsystem_B.Gkpd_c[8]) *
      (Subsystem_B.G_pr_l - Subsystem_B.Gkpd_c[3]) / (Subsystem_B.XB -
      Subsystem_B.Gkpd_c[3]) + Subsystem_B.Gkpd_c[8];
  }

  if (Subsystem_B.G_pr_l < Subsystem_B.Gkpd_c[0]) {
    Subsystem_B.XB = Subsystem_B.Gkpd_c[1];
    if (Subsystem_B.Gkpd_c[0] == Subsystem_B.Gkpd_c[1]) {
      Subsystem_B.XB = Subsystem_B.Gkpd_c[1] + 0.1;
    }

    Subsystem_B.Product = (Subsystem_B.Gkpd_c[6] - Subsystem_B.Gkpd_c[5]) *
      (Subsystem_B.G_pr_l - Subsystem_B.Gkpd_c[0]) / (Subsystem_B.XB -
      Subsystem_B.Gkpd_c[0]) + Subsystem_B.Gkpd_c[5];
  } else {
    guard1 = false;
    if (Subsystem_B.G_pr_l >= Subsystem_B.Gkpd_c[0]) {
      if (Subsystem_B.G_pr_l <= Subsystem_B.Gkpd_c[1]) {
        Subsystem_B.XA = Subsystem_B.Gkpd_c[0];
        Subsystem_B.XB = Subsystem_B.Gkpd_c[1];
        k1 = 5;
        k2 = 6;
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      if (Subsystem_B.XB == Subsystem_B.XA) {
        Subsystem_B.XB += 0.1;
      }

      Subsystem_B.Product = (Subsystem_B.Gkpd_c[k2] - Subsystem_B.Gkpd_c[k1]) *
        (Subsystem_B.G_pr_l - Subsystem_B.XA) / (Subsystem_B.XB - Subsystem_B.XA)
        + Subsystem_B.Gkpd_c[k1];
    }

    guard1 = false;
    if (Subsystem_B.G_pr_l >= Subsystem_B.Gkpd_c[1]) {
      if (Subsystem_B.G_pr_l <= Subsystem_B.Gkpd_c[2]) {
        Subsystem_B.XA = Subsystem_B.Gkpd_c[1];
        Subsystem_B.XB = Subsystem_B.Gkpd_c[2];
        k1 = 6;
        k2 = 7;
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      if (Subsystem_B.XB == Subsystem_B.XA) {
        Subsystem_B.XB += 0.1;
      }

      Subsystem_B.Product = (Subsystem_B.Gkpd_c[k2] - Subsystem_B.Gkpd_c[k1]) *
        (Subsystem_B.G_pr_l - Subsystem_B.XA) / (Subsystem_B.XB - Subsystem_B.XA)
        + Subsystem_B.Gkpd_c[k1];
    }

    guard1 = false;
    if (Subsystem_B.G_pr_l >= Subsystem_B.Gkpd_c[2]) {
      if (Subsystem_B.G_pr_l <= Subsystem_B.Gkpd_c[3]) {
        Subsystem_B.XA = Subsystem_B.Gkpd_c[2];
        Subsystem_B.XB = Subsystem_B.Gkpd_c[3];
        k1 = 7;
        k2 = 8;
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      if (Subsystem_B.XB == Subsystem_B.XA) {
        Subsystem_B.XB += 0.1;
      }

      Subsystem_B.Product = (Subsystem_B.Gkpd_c[k2] - Subsystem_B.Gkpd_c[k1]) *
        (Subsystem_B.G_pr_l - Subsystem_B.XA) / (Subsystem_B.XB - Subsystem_B.XA)
        + Subsystem_B.Gkpd_c[k1];
    }

    guard1 = false;
    if (Subsystem_B.G_pr_l >= Subsystem_B.Gkpd_c[3]) {
      if (Subsystem_B.G_pr_l <= Subsystem_B.Gkpd_c[4]) {
        Subsystem_B.XA = Subsystem_B.Gkpd_c[3];
        Subsystem_B.XB = Subsystem_B.Gkpd_c[4];
        k1 = 8;
        k2 = 9;
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      if (Subsystem_B.XB == Subsystem_B.XA) {
        Subsystem_B.XB += 0.1;
      }

      Subsystem_B.Product = (Subsystem_B.Gkpd_c[k2] - Subsystem_B.Gkpd_c[k1]) *
        (Subsystem_B.G_pr_l - Subsystem_B.XA) / (Subsystem_B.XB - Subsystem_B.XA)
        + Subsystem_B.Gkpd_c[k1];
    }
  }

  /* End of MATLAB Function: '<S7>/��� �����������' */

  /* Fcn: '<S7>/Fcn1' */
  if (Subsystem_B.Divide_o < 0.0) {
    Subsystem_B.Divide_o = -rt_powd_snf(-Subsystem_B.Divide_o,
      0.28571428571428564);
  } else {
    Subsystem_B.Divide_o = rt_powd_snf(Subsystem_B.Divide_o, 0.28571428571428564);
  }

  /* Product: '<S7>/Product2' incorporates:
   *  Fcn: '<S7>/Fcn1'
   *  Fcn: '<S7>/Fcn2'
   *  Product: '<S7>/Divide1'
   */
  Subsystem_B.XA = (1.0 / Subsystem_B.Product * (Subsystem_B.Divide_o - 1.0) +
                    1.0) * Subsystem_B.T_vh;

  /* Step: '<S3>/Step' */
  if (Subsystem_M->Timing.t[0] < Subsystem_P.Step_Time) {
    Subsystem_B.MathFunction = Subsystem_P.Step_Y0;
  } else {
    Subsystem_B.MathFunction = Subsystem_P.Step_YFinal;
  }

  /* End of Step: '<S3>/Step' */

  /* Product: '<S8>/Divide' incorporates:
   *  Fcn: '<S10>/Fcn4'
   *  Gain: '<S7>/Gain2'
   *  Gain: '<S8>/Gain1'
   *  Integrator: '<S12>/Integrator'
   *  Lookup_n-D: '<S8>/1-D Lookup Table'
   *  ManualSwitch: '<S8>/Manual Switch'
   *  Product: '<S10>/Product2'
   *  Product: '<S7>/Product3'
   *  Product: '<S8>/Product'
   *  Sum: '<S10>/Sum'
   *  Sum: '<S7>/Sum'
   *  Sum: '<S8>/Sum'
   */
  Subsystem_B.Divide_i = (((Subsystem_B.Divide_b / (Subsystem_B.Divide_b - 1.0) *
    (Subsystem_X.Integrator_CSTATE - Subsystem_B.Divide) * Subsystem_B.R *
    Subsystem_B.G_pr + Subsystem_B.MathFunction) - Subsystem_P.k_v /
    (Subsystem_P.k_v - 1.0) * Subsystem_P.R_v * (Subsystem_B.XA -
    Subsystem_B.T_vh) * Subsystem_B.G_k) - Subsystem_B.ManualSwitch *
    look1_binlxpw(Subsystem_B.ManualSwitch, Subsystem_P.n_arr,
                  Subsystem_P.k_tr_arr, 1U)) * (1.0 / (0.010966227112321508 *
    Subsystem_P.I * Subsystem_B.ManualSwitch));

  /* Product: '<S12>/Divide' incorporates:
   *  Gain: '<S12>/Gain'
   *  Gain: '<S14>/Gain'
   *  Gain: '<S14>/Gain1'
   *  Integrator: '<S11>/Integrator'
   *  Integrator: '<S12>/Integrator'
   *  Lookup_n-D: '<S13>/i_k'
   *  MinMax: '<S14>/MinMax'
   *  Product: '<S13>/Divide1'
   *  Product: '<S13>/Product'
   *  Product: '<S7>/Product2'
   *  Sum: '<S12>/Sum'
   */
  Subsystem_B.Divide_f = ((fmin(Subsystem_B.Min1, 1.0 / Subsystem_P.L0 *
    Subsystem_B.G_k) * (Subsystem_P.Hu * Subsystem_P.kpd_gor) + look1_binlxpw
    (Subsystem_B.XA, Subsystem_P.T_k_arr, Subsystem_P.i_k_arr, 7U) /
    Subsystem_B.XA * Subsystem_B.G_k * Subsystem_B.XA) - Switch_idx_0) * (1.0 /
    (Subsystem_P.V_ks * Subsystem_X.Integrator_CSTATE_e) * Subsystem_B.R) /
    Divide1_j * Subsystem_X.Integrator_CSTATE;

  /* Sum: '<S11>/Sum1' incorporates:
   *  Gain: '<S11>/Gain'
   *  Integrator: '<S11>/Integrator'
   *  Integrator: '<S12>/Integrator'
   *  Product: '<S11>/Divide'
   *  Product: '<S11>/Divide1'
   *  Sum: '<S11>/Sum'
   */
  Subsystem_B.Sum1 = ((Subsystem_B.Min1 + Subsystem_B.G_k) - Subsystem_B.G_pr) *
    Subsystem_B.R * Subsystem_X.Integrator_CSTATE * (1.0 / Subsystem_P.V_ks) +
    1.0 / Subsystem_X.Integrator_CSTATE * Subsystem_B.Divide_f *
    Subsystem_X.Integrator_CSTATE_e;

  /* Outport: '<Root>/n_pr' incorporates:
   *  Gain: '<S8>/Gain2'
   *  Integrator: '<S8>/Integrator'
   */
  Subsystem_Y.n_pr = Subsystem_P.Gain2_Gain_p * Subsystem_X.Integrator_CSTATE_o;
  if (rtmIsMajorTimeStep(Subsystem_M)) {
    /* Update for Integrator: '<S12>/Integrator' */
    Subsystem_DW.Integrator_IWORK = 0;

    /* Update for Integrator: '<S11>/Integrator' */
    Subsystem_DW.Integrator_IWORK_a = 0;
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(Subsystem_M)) {
    rt_ertODEUpdateContinuousStates(&Subsystem_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++Subsystem_M->Timing.clockTick0;
    Subsystem_M->Timing.t[0] = rtsiGetSolverStopTime(&Subsystem_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      Subsystem_M->Timing.clockTick1++;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void Subsystem_derivatives(void)
{
  XDot_Subsystem_T *_rtXdot;
  _rtXdot = ((XDot_Subsystem_T *) Subsystem_M->derivs);

  /* Derivatives for Integrator: '<S12>/Integrator' */
  _rtXdot->Integrator_CSTATE = Subsystem_B.Divide_f;

  /* Derivatives for Integrator: '<S11>/Integrator' */
  _rtXdot->Integrator_CSTATE_e = Subsystem_B.Sum1;

  /* Derivatives for Integrator: '<S8>/Integrator' */
  _rtXdot->Integrator_CSTATE_o = Subsystem_B.Divide_i;
}

/* Model initialize function */
void Subsystem_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Subsystem_M->solverInfo,
                          &Subsystem_M->Timing.simTimeStep);
    rtsiSetTPtr(&Subsystem_M->solverInfo, &rtmGetTPtr(Subsystem_M));
    rtsiSetStepSizePtr(&Subsystem_M->solverInfo, &Subsystem_M->Timing.stepSize0);
    rtsiSetdXPtr(&Subsystem_M->solverInfo, &Subsystem_M->derivs);
    rtsiSetContStatesPtr(&Subsystem_M->solverInfo, (real_T **)
                         &Subsystem_M->contStates);
    rtsiSetNumContStatesPtr(&Subsystem_M->solverInfo,
      &Subsystem_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&Subsystem_M->solverInfo,
      &Subsystem_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&Subsystem_M->solverInfo,
      &Subsystem_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&Subsystem_M->solverInfo,
      &Subsystem_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&Subsystem_M->solverInfo, (&rtmGetErrorStatus
      (Subsystem_M)));
    rtsiSetRTModelPtr(&Subsystem_M->solverInfo, Subsystem_M);
  }

  rtsiSetSimTimeStep(&Subsystem_M->solverInfo, MAJOR_TIME_STEP);
  Subsystem_M->intgData.x0 = Subsystem_M->odeX0;
  Subsystem_M->intgData.f0 = Subsystem_M->odeF0;
  Subsystem_M->intgData.x1start = Subsystem_M->odeX1START;
  Subsystem_M->intgData.f1 = Subsystem_M->odeF1;
  Subsystem_M->intgData.Delta = Subsystem_M->odeDELTA;
  Subsystem_M->intgData.E = Subsystem_M->odeE;
  Subsystem_M->intgData.fac = Subsystem_M->odeFAC;

  /* initialize */
  {
    int_T i;
    real_T *f = Subsystem_M->intgData.fac;
    for (i = 0; i < (int_T)(sizeof(Subsystem_M->odeFAC)/sizeof(real_T)); i++) {
      f[i] = 1.5e-8;
    }
  }

  Subsystem_M->intgData.DFDX = Subsystem_M->odeDFDX;
  Subsystem_M->intgData.W = Subsystem_M->odeW;
  Subsystem_M->intgData.pivots = Subsystem_M->odePIVOTS;
  Subsystem_M->intgData.xtmp = Subsystem_M->odeXTMP;
  Subsystem_M->intgData.ztmp = Subsystem_M->odeZTMP;
  Subsystem_M->intgData.isFirstStep = true;
  rtsiSetSolverExtrapolationOrder(&Subsystem_M->solverInfo, 4);
  rtsiSetSolverNumberNewtonIterations(&Subsystem_M->solverInfo, 1);
  Subsystem_M->contStates = ((X_Subsystem_T *) &Subsystem_X);
  rtsiSetSolverData(&Subsystem_M->solverInfo, (void *)&Subsystem_M->intgData);
  rtsiSetSolverName(&Subsystem_M->solverInfo,"ode14x");
  rtmSetTPtr(Subsystem_M, &Subsystem_M->Timing.tArray[0]);
  Subsystem_M->Timing.stepSize0 = 0.001;
  rtmSetFirstInitCond(Subsystem_M, 1);

  /* InitializeConditions for Integrator: '<S12>/Integrator' incorporates:
   *  Integrator: '<S11>/Integrator'
   */
  if (rtmIsFirstInitCond(Subsystem_M)) {
    Subsystem_X.Integrator_CSTATE = 0.0;
    Subsystem_X.Integrator_CSTATE_e = 0.0;
  }

  Subsystem_DW.Integrator_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S12>/Integrator' */

  /* InitializeConditions for Integrator: '<S11>/Integrator' */
  Subsystem_DW.Integrator_IWORK_a = 1;

  /* InitializeConditions for Integrator: '<S8>/Integrator' */
  Subsystem_X.Integrator_CSTATE_o = Subsystem_P.Integrator_IC;

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(Subsystem_M)) {
    rtmSetFirstInitCond(Subsystem_M, 0);
  }
}

/* Model terminate function */
void Subsystem_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
