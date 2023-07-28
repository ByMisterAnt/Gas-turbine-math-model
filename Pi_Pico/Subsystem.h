#ifndef RTW_HEADER_Subsystem_h_
#define RTW_HEADER_Subsystem_h_
#include <math.h>
#include <string.h>
#include <stddef.h>
#ifndef Subsystem_COMMON_INCLUDES_
#define Subsystem_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Subsystem_COMMON_INCLUDES_ */

#include "Subsystem_types.h"
#include "rt_nonfinite.h"
#include "math.h"
#include "rt_matrixlib.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeDELTA
#define rtmGetOdeDELTA(rtm)            ((rtm)->odeDELTA)
#endif

#ifndef rtmSetOdeDELTA
#define rtmSetOdeDELTA(rtm, val)       ((rtm)->odeDELTA = (val))
#endif

#ifndef rtmGetOdeDFDX
#define rtmGetOdeDFDX(rtm)             ((rtm)->odeDFDX)
#endif

#ifndef rtmSetOdeDFDX
#define rtmSetOdeDFDX(rtm, val)        ((rtm)->odeDFDX = (val))
#endif

#ifndef rtmGetOdeE
#define rtmGetOdeE(rtm)                ((rtm)->odeE)
#endif

#ifndef rtmSetOdeE
#define rtmSetOdeE(rtm, val)           ((rtm)->odeE = (val))
#endif

#ifndef rtmGetOdeF0
#define rtmGetOdeF0(rtm)               ((rtm)->odeF0)
#endif

#ifndef rtmSetOdeF0
#define rtmSetOdeF0(rtm, val)          ((rtm)->odeF0 = (val))
#endif

#ifndef rtmGetOdeF1
#define rtmGetOdeF1(rtm)               ((rtm)->odeF1)
#endif

#ifndef rtmSetOdeF1
#define rtmSetOdeF1(rtm, val)          ((rtm)->odeF1 = (val))
#endif

#ifndef rtmGetOdeFAC
#define rtmGetOdeFAC(rtm)              ((rtm)->odeFAC)
#endif

#ifndef rtmSetOdeFAC
#define rtmSetOdeFAC(rtm, val)         ((rtm)->odeFAC = (val))
#endif

#ifndef rtmGetOdePIVOTS
#define rtmGetOdePIVOTS(rtm)           ((rtm)->odePIVOTS)
#endif

#ifndef rtmSetOdePIVOTS
#define rtmSetOdePIVOTS(rtm, val)      ((rtm)->odePIVOTS = (val))
#endif

#ifndef rtmGetOdeW
#define rtmGetOdeW(rtm)                ((rtm)->odeW)
#endif

#ifndef rtmSetOdeW
#define rtmSetOdeW(rtm, val)           ((rtm)->odeW = (val))
#endif

#ifndef rtmGetOdeX0
#define rtmGetOdeX0(rtm)               ((rtm)->odeX0)
#endif

#ifndef rtmSetOdeX0
#define rtmSetOdeX0(rtm, val)          ((rtm)->odeX0 = (val))
#endif

#ifndef rtmGetOdeX1START
#define rtmGetOdeX1START(rtm)          ((rtm)->odeX1START)
#endif

#ifndef rtmSetOdeX1START
#define rtmSetOdeX1START(rtm, val)     ((rtm)->odeX1START = (val))
#endif

#ifndef rtmGetOdeXTMP
#define rtmGetOdeXTMP(rtm)             ((rtm)->odeXTMP)
#endif

#ifndef rtmSetOdeXTMP
#define rtmSetOdeXTMP(rtm, val)        ((rtm)->odeXTMP = (val))
#endif

#ifndef rtmGetOdeZTMP
#define rtmGetOdeZTMP(rtm)             ((rtm)->odeZTMP)
#endif

#ifndef rtmSetOdeZTMP
#define rtmSetOdeZTMP(rtm, val)        ((rtm)->odeZTMP = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T Gkpd[16];
  real_T Gkpd_m[14];
  real_T Gkpd_c[10];
  real_T T;                            /* '<S2>/Sum' */
  real_T p;                            /* '<S2>/Gain' */
  real_T Gain;                         /* '<S9>/Gain' */
  real_T T_vh;                         /* '<S5>/Product' */
  real_T Fcn;                          /* '<S7>/Fcn' */
  real_T p_vh;                         /* '<S5>/Gain1' */
  real_T Gain1_m;                      /* '<S7>/Gain1' */
  real_T Divide_i;                     /* '<S8>/Divide' */
  real_T Divide_f;                     /* '<S12>/Divide' */
  real_T Sum1;                         /* '<S11>/Sum1' */
  real_T XA;
  real_T XB;
  real_T MathFunction;                 /* '<S2>/Math Function' */
  real_T R;                            /* '<S6>/R' */
  real_T G_pr;                         /* '<S10>/��������������  �������' */
  real_T G_pr_l;                       /* '<S7>/��������������  �����������' */
  real_T uDLookupTable;                /* '<S8>/1-D Lookup Table' */
  real_T i_k_p;                        /* '<S13>/i_k' */
  real_T kpd;                          /* '<S10>/��� �������' */
  real_T ManualSwitch;
  real_T Divide_b;
  real_T Divide_o;
  real_T G_k;
  real_T Min1;
  real_T Divide;
  real_T Product;
} B_Subsystem_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  int_T Integrator_IWORK;              /* '<S12>/Integrator' */
  int_T Integrator_IWORK_a;            /* '<S11>/Integrator' */
} DW_Subsystem_T;

/* Continuous states (default storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S12>/Integrator' */
  real_T Integrator_CSTATE_e;          /* '<S11>/Integrator' */
  real_T Integrator_CSTATE_o;          /* '<S8>/Integrator' */
} X_Subsystem_T;

/* State derivatives (default storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S12>/Integrator' */
  real_T Integrator_CSTATE_e;          /* '<S11>/Integrator' */
  real_T Integrator_CSTATE_o;          /* '<S8>/Integrator' */
} XDot_Subsystem_T;

/* State disabled  */
typedef struct {
  boolean_T Integrator_CSTATE;         /* '<S12>/Integrator' */
  boolean_T Integrator_CSTATE_e;       /* '<S11>/Integrator' */
  boolean_T Integrator_CSTATE_o;       /* '<S8>/Integrator' */
} XDis_Subsystem_T;

#ifndef ODE14X_INTG
#define ODE14X_INTG

/* ODE14X Integration Data */
typedef struct {
  real_T *x0;
  real_T *f0;
  real_T *x1start;
  real_T *f1;
  real_T *Delta;
  real_T *E;
  real_T *fac;
  real_T *DFDX;
  real_T *W;
  int_T *pivots;
  real_T *xtmp;
  real_T *ztmp;
  real_T *M;
  real_T *M1;
  real_T *Edot;
  real_T *xdot;
  real_T *fminusMxdot;
  boolean_T isFirstStep;
} ODE14X_IntgData;

#endif

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T G_fuel;                       /* '<Root>/G_fuel' */
} ExtU_Subsystem_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T P;                            /* '<Root>/P' */
  real_T n_pr;                         /* '<Root>/n_pr' */
  real_T n;                            /* '<Root>/n' */
} ExtY_Subsystem_T;

/* Parameters (default storage) */
struct P_Subsystem_T_ {
  real_T F_s;                          /* Variable: F_s
                                        * Referenced by: '<S9>/Gain1'
                                        */
  real_T Hu;                           /* Variable: Hu
                                        * Referenced by: '<S14>/Gain1'
                                        */
  real_T I;                            /* Variable: I
                                        * Referenced by: '<S8>/Gain1'
                                        */
  real_T K_arr[12];                    /* Variable: K_arr
                                        * Referenced by: '<S6>/R'
                                        */
  real_T L0;                           /* Variable: L0
                                        * Referenced by: '<S14>/Gain'
                                        */
  real_T R_arr[12];                    /* Variable: R_arr
                                        * Referenced by: '<S6>/R'
                                        */
  real_T R_v;                          /* Variable: R_v
                                        * Referenced by: '<S7>/Gain2'
                                        */
  real_T T_g_arr[8];                   /* Variable: T_g_arr
                                        * Referenced by: '<S15>/i_g'
                                        */
  real_T T_k_arr[8];                   /* Variable: T_k_arr
                                        * Referenced by:
                                        *   '<S13>/i_k'
                                        *   '<S15>/i_k'
                                        */
  real_T V_ks;                         /* Variable: V_ks
                                        * Referenced by:
                                        *   '<S11>/Gain'
                                        *   '<S12>/Gain'
                                        */
  real_T fi;                           /* Variable: fi
                                        * Referenced by: '<S9>/Gain2'
                                        */
  real_T i_g_arr[8];                   /* Variable: i_g_arr
                                        * Referenced by: '<S15>/i_g'
                                        */
  real_T i_k_arr[8];                   /* Variable: i_k_arr
                                        * Referenced by:
                                        *   '<S13>/i_k'
                                        *   '<S15>/i_k'
                                        */
  real_T k_tr_arr[2];                  /* Variable: k_tr_arr
                                        * Referenced by: '<S8>/1-D Lookup Table'
                                        */
  real_T k_v;                          /* Variable: k_v
                                        * Referenced by: '<S7>/Gain2'
                                        */
  real_T kpd_gor;                      /* Variable: kpd_gor
                                        * Referenced by: '<S14>/Gain1'
                                        */
  real_T n_arr[2];                     /* Variable: n_arr
                                        * Referenced by: '<S8>/1-D Lookup Table'
                                        */
  real_T sig_g;                        /* Variable: sig_g
                                        * Referenced by: '<S11>/Gain1'
                                        */
  real_T sig_s;                        /* Variable: sig_s
                                        * Referenced by: '<S9>/Gain'
                                        */
  real_T sig_vh;                       /* Variable: sig_vh
                                        * Referenced by: '<S5>/Gain1'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 90500
                                        * Referenced by: '<S8>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: 1
                                        * Referenced by: '<S8>/Saturation'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: 10000000
                                        * Referenced by: '<S8>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: 1
                                        * Referenced by: '<S8>/Saturation1'
                                        */
  real_T n_max_Value;                  /* Expression: 0.0084
                                        * Referenced by: '<S4>/n_max'
                                        */
  real_T Constant_Value;               /* Expression: 0.002
                                        * Referenced by: '<S4>/Constant'
                                        */
  real_T Constant_Value_g;             /* Expression: 288
                                        * Referenced by: '<S2>/Constant'
                                        */
  real_T _Value;                       /* Expression: 1000
                                        * Referenced by: '<S1>/������ (�)'
                                        */
  real_T Gain1_Gain;                   /* Expression: 6.5/1000
                                        * Referenced by: '<S2>/Gain1'
                                        */
  real_T Gain2_Gain;                   /* Expression: 6.5/1000/288
                                        * Referenced by: '<S2>/Gain2'
                                        */
  real_T Constant1_Value;              /* Expression: 1
                                        * Referenced by: '<S2>/Constant1'
                                        */
  real_T Constant2_Value;              /* Expression: 5.25
                                        * Referenced by: '<S2>/Constant2'
                                        */
  real_T Gain_Gain;                    /* Expression: 1.013e5
                                        * Referenced by: '<S2>/Gain'
                                        */
  real_T Gain2_Gain_j;                 /* Expression: 1/101325
                                        * Referenced by: '<S10>/Gain2'
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<S8>/Integrator'
                                        */
  real_T Gain_Gain_h;                  /* Expression: 1/60
                                        * Referenced by: '<S8>/Gain'
                                        */
  real_T Switch_Threshold;             /* Expression: 0
                                        * Referenced by: '<S15>/Switch'
                                        */
  real_T _Value_a;                     /* Expression: 0
                                        * Referenced by: '<S1>/�������� (���)'
                                        */
  real_T Saturation_UpperSat_k;        /* Expression: 0.8
                                        * Referenced by: '<S7>/Saturation'
                                        */
  real_T Saturation_LowerSat_b;        /* Expression: 0.05
                                        * Referenced by: '<S7>/Saturation'
                                        */
  real_T Gain1_Gain_b;                 /* Expression: 1/101325
                                        * Referenced by: '<S7>/Gain1'
                                        */
  real_T Step_Time;                    /* Expression: 3
                                        * Referenced by: '<S3>/Step'
                                        */
  real_T Step_Y0;                      /* Expression: 60000
                                        * Referenced by: '<S3>/Step'
                                        */
  real_T Step_YFinal;                  /* Expression: 0
                                        * Referenced by: '<S3>/Step'
                                        */
  real_T Gain2_Gain_p;                 /* Expression: 1/60
                                        * Referenced by: '<S8>/Gain2'
                                        */
  uint8_T ManualSwitch_CurrentSetting;
                              /* Computed Parameter: ManualSwitch_CurrentSetting
                               * Referenced by: '<S8>/Manual Switch'
                               */
};

/* Real-time Model Data Structure */
struct tag_RTM_Subsystem_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_Subsystem_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeX0[3];
  real_T odeF0[3];
  real_T odeX1START[3];
  real_T odeF1[3];
  real_T odeDELTA[3];
  real_T odeE[4*3];
  real_T odeFAC[3];
  real_T odeDFDX[3*3];
  real_T odeW[3*3];
  int_T odePIVOTS[3];
  real_T odeXTMP[3];
  real_T odeZTMP[3];
  ODE14X_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    boolean_T firstInitCondFlag;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
extern P_Subsystem_T Subsystem_P;

/* Block signals (default storage) */
extern B_Subsystem_T Subsystem_B;

/* Continuous states (default storage) */
extern X_Subsystem_T Subsystem_X;

/* Block states (default storage) */
extern DW_Subsystem_T Subsystem_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_Subsystem_T Subsystem_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_Subsystem_T Subsystem_Y;

/* Model entry point functions */
extern void Subsystem_initialize(void);
extern void Subsystem_step(void);
extern void Subsystem_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Subsystem_T *const Subsystem_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S7>/Scope3' : Unused code path elimination
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Note that this particular code originates from a subsystem build,
 * and has its own system numbers different from the parent model.
 * Refer to the system hierarchy for this subsystem below, and use the
 * MATLAB hilite_system command to trace the generated code back
 * to the parent model.  For example,
 *
 * hilite_system('dvig_16_true/Subsystem')    - opens subsystem dvig_16_true/Subsystem
 * hilite_system('dvig_16_true/Subsystem/Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'dvig_16_true'
 * '<S1>'   : 'dvig_16_true/Subsystem'
 * '<S2>'   : 'dvig_16_true/Subsystem/���������'
 * '<S3>'   : 'dvig_16_true/Subsystem/����'
 * '<S4>'   : 'dvig_16_true/Subsystem/��������� ������� �������'
 * '<S5>'   : 'dvig_16_true/Subsystem/����/������� ����������'
 * '<S6>'   : 'dvig_16_true/Subsystem/����/������ ��������'
 * '<S7>'   : 'dvig_16_true/Subsystem/����/����������'
 * '<S8>'   : 'dvig_16_true/Subsystem/����/�����'
 * '<S9>'   : 'dvig_16_true/Subsystem/����/�����'
 * '<S10>'  : 'dvig_16_true/Subsystem/����/�������'
 * '<S11>'  : 'dvig_16_true/Subsystem/����/������ ��������/��������'
 * '<S12>'  : 'dvig_16_true/Subsystem/����/������ ��������/�����������'
 * '<S13>'  : 'dvig_16_true/Subsystem/����/������ ��������/�����������/������� ����. �� �����'
 * '<S14>'  : 'dvig_16_true/Subsystem/����/������ ��������/�����������/������� ��������'
 * '<S15>'  : 'dvig_16_true/Subsystem/����/������ ��������/�����������/���� ������� ��'
 * '<S16>'  : 'dvig_16_true/Subsystem/����/����������/��� �����������'
 * '<S17>'  : 'dvig_16_true/Subsystem/����/����������/��������������  �����������'
 * '<S18>'  : 'dvig_16_true/Subsystem/����/�������/��� �������'
 * '<S19>'  : 'dvig_16_true/Subsystem/����/�������/��������������  �������'
 */
#endif                                 /* RTW_HEADER_Subsystem_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
