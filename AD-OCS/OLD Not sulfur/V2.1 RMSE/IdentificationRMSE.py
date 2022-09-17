# IDENTIFICATION


from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt

from minimRMSE import minimopt

from dataimport import*
from PhysConstants import*
from ReactorConf import*

Dil     = 1/HRT

# KINETICS

Xr1 = [Dil*S1, Dil]
Yr1 = S1

def func1(X,a,b):
    x1,x2 = X
    return a*x1 +a*b*x2 + b*C_d[0]/(1-C_d[0])

popt1, pcov1 = curve_fit(func1, Xr1, Yr1)
KS1 = popt1[1]                          # [1/d]      Max biomass growth rate
mu1_max = alfa/(1-C_d[0])/popt1[0]      # [g/L]      Half saturation constant

Xr2 = (Dil, S2)
Yr2 = S2

def func2(X,a,b,c):
    x1,x2 = X
    return a*x1*x2 + a*b*x1 + C_d[1]/(1-C_d[1])*b + a*c*x1*(x2**2) + C_d[1]/(1-C_d[1])*c*(x2**2)

guess2 = [10, 3, 1/300]

popt2, pcov2 = curve_fit(func2, Xr2, Yr2, guess2)
# res_lsq2 = least_squares(fun, x0, args=(t_train, y_train))

KS2     = popt2[1];
mu2_max = alfa/(1-C_d[1])/popt2[0];
KI2     = 1/popt2[2];


mu_max = (mu1_max, mu2_max)         # [1/d]      Max biomass growth rate
Ks     = (KS1, KS2)                 # [g/L]      Half saturation constant
KI2    = np.abs(KI2)                # [mmol/L]   Inhibition constant for S2
kd     = np.multiply(C_d,mu_max)

# Hydrolysis

X_hydr = XT
Y_hydr = Dil*(XT_in-XT)
res_hydr = linregress(X_hydr, Y_hydr)
k_hyd = res_hydr.slope

# Physical - L/G TRANSFER
pKb = -np.log10(Kb)
fun = 1/(1+10**(pH-pKb))
X3r = C*fun - KH*P_C
Y3r = q_C

mdl3 = linregress(X3r, Y3r)
kLa = mdl3.slope                      # [1/d]      L/G transfer rate

# Yield coefficients: [CODdeg, VFAprof, VFAcons, CO2prod(1), CO2prod(2), CH4prod, hydrolysis]
mu_1 = alfa*Dil + C_d[0]*mu1_max
mu_2 = alfa*Dil + C_d[1]*mu2_max


minsol_4 = minimopt(4, k_hyd, mu_1, mu_2)
k1   = float(minsol_4)

minsol_5 = minimopt(5, k_hyd, mu_1, mu_2)
k6   = float(minsol_5)

# DOUBLE REGRESSIONS

minsol_6 = minimopt(6, k_hyd, mu_1, mu_2)
AA  = minsol_6[0]
BB  = minsol_6[1]

minsol_7 = minimopt(7, k_hyd,mu_1, mu_2)
CC = minsol_7[0]
DD = minsol_7[1]

k2 = BB*k1
k3 = k6/AA
k4 = CC*k1
k5 = DD*k6

k = [k1, k2, k3, k4, k5, k6, k_hyd]
