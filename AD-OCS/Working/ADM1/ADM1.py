import math
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

f_sI_xc = 0.1
f_xI_xc = 0.2
f_ch_xc = 0.2
f_pr_xc = 0.2
f_li_xc = 0.3
N_xc = 0.002685714
N_I = 0.004285714
N_aa = 0.007
C_xc = 0.02786
C_sI = 0.03
C_ch = 0.0313
C_pr = 0.03
C_li = 0.022
C_xI = 0.03
C_su = 0.0313
C_aa = 0.03
f_fa_li = 0.95
C_fa = 0.0217
f_h2_su = 0.19
f_bu_su = 0.13
f_pro_su = 0.27
f_ac_su = 0.41
N_bac = 0.005714286
C_bu = 0.025
C_pro = 0.0268
C_ac = 0.0313
C_bac = 0.0313
Y_su = 0.1
f_h2_aa = 0.06
f_va_aa = 0.23
f_bu_aa = 0.26
f_pro_aa = 0.05
f_ac_aa = 0.4
C_va = 0.024
Y_aa = 0.08
Y_fa = 0.06
Y_c4 = 0.06
Y_pro = 0.04
C_ch4 = 0.0156
Y_ac = 0.05
Y_h2 = 0.06
f_ac_li = 0.7
f_h2_li = 0.3
f_pro_va = 0.54
f_ac_va = 0.31
f_h2_va = 0.15
f_ac_bu = 0.8
f_h2_bu = 0.2
f_ac_pro = 0.57
f_h2_pro = 0.43

k_dis = 0.5
k_hyd_ch = 10
k_hyd_pr = 10
k_hyd_li = 10
K_S_IN = 0.0001
km_su = 30
ks_su = 0.5
pH_UL_aa = 5.5
pH_LL_aa = 4
km_aa = 50
ks_aa = 0.3
km_fa = 6
ks_fa = 0.4
K_Ih2_fa = 0.000005
km_c4 = 20
ks_c4 = 0.2
K_Ih2_c4 = 0.00001
km_pro = 13
ks_pro = 0.1
K_Ih2_pro = 0.0000035
km_ac = 8
ks_ac = 0.15
K_I_nh3 = 0.0018
pH_UL_ac = 7
pH_LL_ac = 6
km_h2 = 35
ks_h2 = 0.000007
pH_UL_h2 = 6
pH_LL_h2 = 5
k_dec_Xsu = 0.02
k_dec_Xaa = 0.02
k_dec_Xfa = 0.02
k_dec_Xc4 = 0.02
k_dec_Xpro = 0.02
k_dec_Xac = 0.02
k_dec_Xh2 = 0.02

V_liq = 3400 # ========================================================================== Volumes
V_gas = 300

R_gas = 0.083145  # universal gas constant dm3*bar/(mol*K) = 8.3145 J/(mol*K)
T_base = 298.15
T_op = 308.15 # ========================================================================== Temperature
K_w = math.exp(55900/(R_gas*100)*(1/T_base-1/T_op))*10**(-14);
K_a_va = 1.38*10**(-5)
K_a_bu = 1.51*10**(-5)
K_a_pro = 1.32*10**(-5)
K_a_ac = 1.74*10**(-5)
K_a_co2 = 10**(-6.35)*math.exp(7646/(R_gas*100)*(1/T_base-1/T_op))
K_a_IN = 10**(-9.25)*math.exp(51965/(R_gas*100)*(1/T_base-1/T_op))
k_A_Bva = 1e10
k_A_Bbu = 1e10
k_A_Bpro = 1e10
k_A_Bac = 1e10
k_A_B_co2 = 1e10
k_A_BIN = 1e10
p_atm = 1.01325
p_gas_h2o = 0.0313*np.exp(5290*(1/T_base-1/T_op))
k_p = 50000
k_L_a = 200
K_H_co2 = 0.035*math.exp(-19140/(R_gas*100)*(1/T_base-1/T_op))
K_H_ch4 = 0.0014*math.exp(-14240/(R_gas*100)*(1/T_base-1/T_op))
K_H_h2 = 0.00078*math.exp(-4180/(R_gas*100)*(1/T_base-1/T_op))

W_xc_in =  2
W_ch_in =  5
W_pr_in =  20
W_li_in =  5
S_su_in = 0.01 #
X_su_in = 0     
S_aa_in = 0.001 # 
X_aa_in = 0.01 
S_fa_in = 0.001 #
X_fa_in = 0.01
S_va_in = 0.001 #
S_bu_in = 0.001
X_c4_in = 0.01
S_pro_in = 0.001
X_pro_in = 0.01
S_ac_in = 0.001
X_ac_in = 0.01
S_h2_in = 10**(-8)
X_h2_in = 0.01
S_ch4_in = 10**(-5)
S_va_dis_in = 0
S_bu_dis_in = 0
S_pro_dis_in = 0
S_ac_dis_in = 0
S_hco3_in = 0
S_IC_in = 0.04
S_gas_h2_in = 0
p_gas_h2_in = 0
S_gas_co2_in = 0
p_gas_co2_in = 0
S_gas_ch4_in = 0
p_gas_ch4_in = 0
S_Hplus_in = 0
S_IN_in = 0.01
S_nh3_in = 0
S_nh4_in = 0
q_in = 170
S_I_in = 0.02
W_I_in = 25
S_cat_in = 0.04
S_an_in = 0.02

def adm1_ss(x,t):

    W_xc = x[0]
    W_ch = x[1]
    W_pr = x[2]
    W_li = x[3]
    S_su = x[4]
    X_su = x[5]
    S_aa = x[6]
    X_aa = x[7]
    S_fa = x[8]
    X_fa = x[9]
    S_va = x[10]
    S_bu = x[11]
    X_c4 = x[12]
    S_pro = x[13]
    X_pro = x[14]
    S_ac = x[15]
    X_ac = x[16]
    S_h2 = x[17]
    X_h2 = x[18]
    S_ch4 = x[19]
    S_va_dis = x[20]
    S_bu_dis = x[21]
    S_pro_dis = x[22]
    S_ac_dis = x[23]
    S_hco3 = x[24]
    S_IC = x[25]
    S_gas_h2 = x[26]
    S_gas_ch4 = x[27]
    S_gas_co2 = x[28]
    S_IN = x[29]
    S_nh3 = x[30]
    S_I = x[31]
    W_I = x[32]
    S_cat = x[33]
    S_an = x[34]

    S_co2 = S_IC-S_hco3
    S_nh4 = S_IN-S_nh3
    theta = S_cat+S_nh4-S_hco3-S_ac_dis/64-S_pro_dis/112-S_bu_dis/160-S_va_dis/208-S_an

    S_Hplus= (-theta/2+0.5*np.sqrt(theta**2+4*K_w))

    if S_Hplus > 0.0001:
        S_Hplus = 0.00001

    p_H1 = -math.log(S_Hplus,10)
    pHstd = 7;

    if p_H1 < 5:
        p_H = pHstd;
    else:
        p_H = p_H1;

    KpH_aa = 10**((-pH_LL_aa+pH_UL_aa)/2)
    KpH_ac = 10**((-pH_LL_ac+pH_UL_ac)/2)
    KpH_h2 = 10**((-pH_LL_h2+pH_UL_h2)/2)

    I_pH_aa = KpH_aa**2/(S_Hplus**2+KpH_aa**2)
    I_pH_ac = KpH_ac**3/(S_Hplus**3+KpH_ac**3)
    I_pH_h2 = KpH_h2**3/(S_Hplus**3+KpH_h2**3)

    I_IN_lim = 1/(1+K_S_IN/S_IN)
    I_h2_fa = 1/(1+S_h2/K_Ih2_fa)
    I_h2_c4 = 1/(1+S_h2/K_Ih2_c4)
    I_h2_pro = 1/(1+S_h2/K_Ih2_pro)
    I_nh3 = 1/(1+S_nh3/K_I_nh3)
    I5 = I_pH_aa*I_IN_lim
    I6 = I_pH_aa*I_IN_lim
    I7 = I_pH_aa*I_IN_lim*I_h2_fa
    I8 = I_pH_aa*I_IN_lim*I_h2_c4
    I9 = I_pH_aa*I_IN_lim*I_h2_c4
    I10 = I_pH_aa*I_IN_lim*I_h2_pro
    I11 = I_pH_ac*I_IN_lim*I_nh3
    I12 = I_pH_h2*I_IN_lim

    rho1 = k_dis*W_xc
    rho2 = k_hyd_ch*W_ch
    rho3 = k_hyd_pr*W_pr
    rho4 = k_hyd_li*W_li
    rho5 = km_su*S_su/(ks_su+S_su)*X_su*I5
    rho6 = km_aa*S_aa/(ks_aa+S_aa)*X_aa*I6
    rho7 = km_fa*S_fa/(ks_fa+S_fa)*X_fa*I7
    rho8 = km_c4*S_va/(ks_c4+S_va)*S_va/(S_bu+S_va+10**(-6))*X_c4*I8
    rho9 = km_c4*S_bu/(ks_c4+S_bu)*S_bu/(S_bu+S_va+10**(-6))*X_c4*I9
    rho10 = km_pro*S_pro/(ks_pro+S_pro)*X_pro*I10
    rho11 = km_ac*S_ac/(ks_ac+S_ac)*X_ac*I11
    rho12 = km_h2*S_h2/(ks_h2+S_h2)*X_h2*I12
    rho13 = k_dec_Xsu*X_su
    rho14 = k_dec_Xaa*X_aa
    rho15 = k_dec_Xfa*X_fa
    rho16 = k_dec_Xc4*X_c4
    rho17 = k_dec_Xpro*X_pro
    rho18 = k_dec_Xac*X_ac
    rho19 = k_dec_Xh2*X_h2

    p_gas_h2 = S_gas_h2*R_gas*T_op/16
    p_gas_ch4 = S_gas_ch4*R_gas*T_op/64
    p_gas_co2 = S_gas_co2*R_gas*T_op
    p_gas = p_gas_h2+p_gas_ch4+p_gas_co2+p_gas_h2o
    q_gas = k_p*(p_gas-p_atm)*p_gas/p_atm

    s_1 = -C_xc+f_sI_xc*C_sI+f_ch_xc*C_ch+f_pr_xc*C_pr+f_li_xc*C_li+f_xI_xc*C_xI
    s_2 = -C_ch+C_su
    s_3 = -C_pr+C_aa
    s_4 = -C_li+(1-f_fa_li)*C_su+f_fa_li*C_fa
    s_5 = -C_su+(1-Y_su)*(f_bu_su*C_bu+f_pro_su*C_pro+f_ac_su*C_ac)+Y_su*C_bac
    s_6 = -C_aa+(1-Y_aa)*(f_va_aa*C_va+f_bu_aa*C_bu+f_pro_aa*C_pro+f_ac_aa*C_ac)+Y_aa*C_bac
    s_7 = -C_fa+(1-Y_fa)*0.7*C_ac+Y_fa*C_bac
    s_8 = -C_va+(1-Y_c4)*0.54*C_pro+(1-Y_c4)*0.31*C_ac+Y_c4*C_bac
    s_9 = -C_bu+(1-Y_c4)*0.8*C_ac+Y_c4*C_bac
    s_10 = -C_pro+(1-Y_pro)*0.57*C_ac+Y_pro*C_bac
    s_11 = -C_ac+(1-Y_ac)*C_ch4 +Y_ac*C_bac
    s_12 = (1-Y_h2)*C_ch4+Y_h2*C_bac
    s_13 = -C_bac+C_xc

    rhoA4 = k_A_Bva*(S_va_dis*(K_a_va+S_Hplus)-K_a_va*S_va)
    rhoA5 = k_A_Bbu*(S_bu_dis*(K_a_bu+S_Hplus)-K_a_bu*S_bu)
    rhoA6 = k_A_Bpro*(S_pro_dis*(K_a_pro+S_Hplus)-K_a_pro*S_pro)
    rhoA7 = k_A_Bac*(S_ac_dis*(K_a_ac+S_Hplus)-K_a_ac*S_ac)
    rhoA10 = k_A_B_co2*(S_hco3*(K_a_co2+S_Hplus)-K_a_co2*S_IC)
    rhoA11 = k_A_BIN*(S_nh3*(K_a_IN+S_Hplus)-K_a_IN*S_IN)

    rhoT8 = k_L_a*(S_h2-16*K_H_h2*p_gas_h2)
    rhoT9 = k_L_a*(S_ch4-64*K_H_ch4*p_gas_ch4)
    rhoT10 = k_L_a*(S_co2-K_H_co2*p_gas_co2)


    y = np.empty(35)

    y[0] = q_in/V_liq*(W_xc_in-W_xc)-rho1+rho13+rho14+rho15+rho16+rho17+rho18+rho19

    y[1] = q_in/V_liq*(W_ch_in-W_ch)+f_ch_xc*rho1-rho2

    y[2] = q_in/V_liq*(W_pr_in-W_pr)+f_pr_xc*rho1-rho3

    y[3] = q_in/V_liq*(W_li_in-W_li)+f_li_xc*rho1-rho4

    y[4] = q_in/V_liq*(S_su_in-S_su)+rho2+(1-f_fa_li)*rho4-rho5

    y[5] = q_in/V_liq*(X_su_in-X_su)+Y_su*rho5-rho13

    y[6] = q_in/V_liq*(S_aa_in-S_aa)+rho3-rho6

    y[7] = q_in/V_liq*(X_aa_in-X_aa)+Y_aa*rho6-rho14

    y[8] = q_in/V_liq*(S_fa_in-S_fa)+f_fa_li*rho4-rho7

    y[9] = q_in/V_liq*(X_fa_in-X_fa)+Y_fa*rho7-rho15

    y[10] = q_in/V_liq*(S_va_in-S_va)+(1-Y_aa)*f_va_aa*rho6-rho8

    y[11] = q_in/V_liq*(S_bu_in-S_bu)+(1-Y_su)*f_bu_su*rho5+(1-Y_aa)*f_bu_aa*rho6-rho9

    y[12] = q_in/V_liq*(X_c4_in-X_c4)+Y_c4*rho8+Y_c4*rho9- rho16

    y[13] = q_in/V_liq*(S_pro_in-S_pro)+(1-Y_su)*f_pro_su*rho5+(1-Y_aa)*f_pro_aa*rho6+(1-Y_c4)*f_pro_va*rho8-rho10

    y[14] = q_in/V_liq*(X_pro_in-X_pro)+Y_pro*rho10-rho17

    y[15] = q_in/V_liq*(S_ac_in-S_ac)+(1-Y_su)*f_ac_su*rho5+(1-Y_aa)*f_ac_aa*rho6+(1-Y_fa)*f_ac_li*rho7+(1-Y_c4)*f_ac_va*rho8+(1-Y_c4)*f_ac_bu*rho9+(1-Y_pro)*f_ac_pro*rho10-rho11

    y[16] = q_in/V_liq*(X_ac_in-X_ac)+Y_ac*rho11-rho18

    y[17] = q_in/V_liq*(S_h2_in-S_h2)+(1-Y_su)*f_h2_su*rho5+(1-Y_aa)*f_h2_aa*rho6+(1-Y_fa)*f_h2_li*rho7+(1-Y_c4)*f_h2_va*rho8+(1-Y_c4)*f_h2_bu*rho9+(1-Y_pro)*f_h2_pro*rho10-rho12-rhoT8

    y[18] = q_in/V_liq*(X_h2_in-X_h2)+Y_h2*rho12-rho19

    y[19] = q_in/V_liq*(S_ch4_in-S_ch4)+(1-Y_ac)*rho11+(1-Y_h2)*rho12-rhoT9

    y[20] = -rhoA4

    y[21] = -rhoA5

    y[22] = -rhoA6

    y[23] = -rhoA7

    y[24] = -rhoA10

    y[25] = q_in/V_liq*(S_IC_in-S_IC)-(s_1*rho1+s_2*rho2+s_3*rho3+s_4*rho4+s_5*rho5+s_6*rho6+s_7*rho7+s_8*rho8+s_9*rho9+s_10*rho10+s_11*rho11+s_12*rho12+s_13*(rho13+rho14+rho15+rho16+rho17+rho18+rho19))-rhoT10

    y[26] = -S_gas_h2*q_gas/V_gas+rhoT8*V_liq/V_gas

    y[27] = -S_gas_ch4*q_gas/V_gas+rhoT9*V_liq/V_gas

    y[28] = -S_gas_co2*q_gas/V_gas+rhoT10*V_liq/V_gas

    y[29] = q_in/V_liq*(S_IN_in-S_IN)-Y_su*N_bac*rho5+(N_aa-Y_aa*N_bac)*rho6-Y_fa*N_bac*rho7-Y_c4*N_bac*rho8-Y_c4*N_bac*rho9-Y_pro*N_bac*rho10-Y_ac*N_bac*rho11-Y_h2*N_bac*rho12+(N_bac-N_xc)*(rho13+rho14+rho15+rho16+rho17+rho18+rho19)+(N_xc-f_xI_xc*N_I-f_sI_xc*N_I-f_pr_xc*N_aa)*rho1

    y[30] = -rhoA11

    y[31] = q_in/V_liq*(S_I_in-S_I)+f_sI_xc*rho1

    y[32] = q_in/V_liq*(W_I_in-W_I)+f_xI_xc*rho1

    y[33] = q_in/V_liq*(S_cat_in-S_cat)

    y[34] = q_in/V_liq*(S_an_in-S_an)

    return y


x0 = np.array( [2, 5, 20, 5, 0.01, 0 ,0.001, 0.01 ,0.001, 0.01, 0.001, 0.001, 0.01 ,0.001, 0.01 ,0.001, 0.01 ,1e-8, 0.01, 1e-5 , 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0.01, 0, 0.02, 25, 0.04, 0.02])
tspan = np.linspace(0,100,100000)

x = odeint(adm1_ss, x0, tspan,atol = 1e-7, rtol = 1e-8, mxstep=10000)

x_df = pd.DataFrame(x[-1,:])

S_co2 = x[:,25]-x[:,24]
Biogas = S_co2+x[:,19]+x[:,17]

S_gas_H2 = x[:,26]
S_gas_Ch4 = x[:,27]
S_gas_Co2 = x[:,28]

P_gas_h2 = S_gas_H2*R_gas*T_op/16
P_gas_ch4 = S_gas_Ch4*R_gas*T_op/64
P_gas_co2 = S_gas_Co2*R_gas*T_op
P_gas = P_gas_h2+P_gas_ch4+P_gas_co2+p_gas_h2o

