#include "Subsystem.h"
#include "Subsystem_private.h"

/* Block parameters (default storage) */
P_Subsystem_T Subsystem_P = {
  /* Variable: F_s
   * Referenced by: '<S9>/Gain1'
   */
  0.004,

  /* Variable: Hu
   * Referenced by: '<S14>/Gain1'
   */
  4.2915E+7,

  /* Variable: I
   * Referenced by: '<S8>/Gain1'
   */
  0.00088,

  /* Variable: K_arr
   * Referenced by: '<S6>/R'
   */
  { 0.0, 1.48, 2.96, 4.44, 5.92, 7.4, 8.88, 10.36, 11.84, 13.32, 14.8, 22.2 },

  /* Variable: L0
   * Referenced by: '<S14>/Gain'
   */
  14.627,

  /* Variable: R_arr
   * Referenced by: '<S6>/R'
   */
  { 510.7, 463.9, 424.2, 390.9, 363.4, 341.2, 323.5, 310.9, 301.3, 294.0, 289.6,
    288.2 },

  /* Variable: R_v
   * Referenced by: '<S7>/Gain2'
   */
  287.0,

  /* Variable: T_g_arr
   * Referenced by: '<S15>/i_g'
   */
  { 373.0, 473.0, 573.0, 673.0, 773.0, 873.0, 973.0, 1073.0 },

  /* Variable: T_k_arr
   * Referenced by:
   *   '<S13>/i_k'
   *   '<S15>/i_k'
   */
  { 223.15, 273.15, 323.15, 373.15, 423.15, 473.15, 523.15, 573.15 },

  /* Variable: V_ks
   * Referenced by:
   *   '<S11>/Gain'
   *   '<S12>/Gain'
   */
  0.0017671458676442591,

  /* Variable: fi
   * Referenced by: '<S9>/Gain2'
   */
  0.98,

  /* Variable: i_g_arr
   * Referenced by: '<S15>/i_g'
   */
  { 377410.0, 480990.0, 586620.0, 694690.0, 805340.0, 918560.0,
    1.0340999999999999E+6, 1.1518E+6 },

  /* Variable: i_k_arr
   * Referenced by:
   *   '<S13>/i_k'
   *   '<S15>/i_k'
   */
  { 222970.0, 273050.0, 323230.0, 373550.0, 424140.0, 475090.0, 524390.0,
    578430.0 },

  /* Variable: k_tr_arr
   * Referenced by: '<S8>/1-D Lookup Table'
   */
  { 0.0174, 0.087 },

  /* Variable: k_v
   * Referenced by: '<S7>/Gain2'
   */
  1.4,

  /* Variable: kpd_gor
   * Referenced by: '<S14>/Gain1'
   */
  0.98,

  /* Variable: n_arr
   * Referenced by: '<S8>/1-D Lookup Table'
   */
  { 0.0, 90457.0 },

  /* Variable: sig_g
   * Referenced by: '<S11>/Gain1'
   */
  0.955,

  /* Variable: sig_s
   * Referenced by: '<S9>/Gain'
   */
  0.95,

  /* Variable: sig_vh
   * Referenced by: '<S5>/Gain1'
   */
  0.99,

  /* Expression: 90500
   * Referenced by: '<S8>/Saturation'
   */
  90500.0,

  /* Expression: 1
   * Referenced by: '<S8>/Saturation'
   */
  1.0,

  /* Expression: 10000000
   * Referenced by: '<S8>/Saturation1'
   */
  1.0E+7,

  /* Expression: 1
   * Referenced by: '<S8>/Saturation1'
   */
  1.0,

  /* Expression: 0.0084
   * Referenced by: '<S4>/n_max'
   */
  0.0084,

  /* Expression: 0.002
   * Referenced by: '<S4>/Constant'
   */
  0.002,

  /* Expression: 288
   * Referenced by: '<S2>/Constant'
   */
  288.0,

  /* Expression: 1000
   * Referenced by: '<S1>/������ (�)'
   */
  1000.0,

  /* Expression: 6.5/1000
   * Referenced by: '<S2>/Gain1'
   */
  0.0065,

  /* Expression: 6.5/1000/288
   * Referenced by: '<S2>/Gain2'
   */
  2.2569444444444443E-5,

  /* Expression: 1
   * Referenced by: '<S2>/Constant1'
   */
  1.0,

  /* Expression: 5.25
   * Referenced by: '<S2>/Constant2'
   */
  5.25,

  /* Expression: 1.013e5
   * Referenced by: '<S2>/Gain'
   */
  101300.0,

  /* Expression: 1/101325
   * Referenced by: '<S10>/Gain2'
   */
  9.8692326671601285E-6,

  /* Expression: 0
   * Referenced by: '<S8>/Integrator'
   */
  0.0,

  /* Expression: 1/60
   * Referenced by: '<S8>/Gain'
   */
  0.016666666666666666,

  /* Expression: 0
   * Referenced by: '<S15>/Switch'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S1>/�������� (���)'
   */
  0.0,

  /* Expression: 0.8
   * Referenced by: '<S7>/Saturation'
   */
  0.8,

  /* Expression: 0.05
   * Referenced by: '<S7>/Saturation'
   */
  0.05,

  /* Expression: 1/101325
   * Referenced by: '<S7>/Gain1'
   */
  9.8692326671601285E-6,

  /* Expression: 3
   * Referenced by: '<S3>/Step'
   */
  3.0,

  /* Expression: 60000
   * Referenced by: '<S3>/Step'
   */
  60000.0,

  /* Expression: 0
   * Referenced by: '<S3>/Step'
   */
  0.0,

  /* Expression: 1/60
   * Referenced by: '<S8>/Gain2'
   */
  0.016666666666666666,

  /* Computed Parameter: ManualSwitch_CurrentSetting
   * Referenced by: '<S8>/Manual Switch'
   */
  0U
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
