clear all
clc

k_v = 1.4;

R_v = 287;
I = 8.8e-4;
l_ks = 0.1;
d_vn = 0.05;
d_nar = 0.2;
F_s = 0.004;
sig_vh = 0.99;
sig_g = 0.955;
sig_s = 0.95;
fi = 0.98;
sig_m = 0.98;
kpd_gor = 0.98;
Hu = 42.915e6;
L0 = 14.627;

%init
F_ks = pi*(d_nar - d_vn)^2/4;
V_ks = F_ks*l_ks;

%������� ��� 1-D ������������
K_arr = [0 1.48 2.96 4.44 5.92 7.40 8.88 10.36 11.84 13.32 14.80 22.20];
R_arr = [510.7 463.9 424.2 390.9 363.4 341.2 323.5 310.9 301.3 294.0 289.6 288.2];

T_k_arr = [223.15 273.15 323.15 373.15 423.15 473.15 523.15 573.15];
i_k_arr = [222.97 273.05 323.23 373.55 424.14 475.09 524.39 578.43].*1000;

T_g_arr = [373 473 573 673 773 873 973 1073];
i_g_arr =[377.41 480.99 586.62 694.69 805.34 918.56 1034.1 1151.8].*1000;

n_arr = [0 90457];
k_tr_arr = [0.0174 0.087];

