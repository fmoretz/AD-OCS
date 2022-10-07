# IDENTIFICATION
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

import pylops
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import ( LinearRegression,
    TheilSenRegressor,
    RANSACRegressor,
    HuberRegressor,
    )

from dataimport import*
from PhysConstants import*
from ReactorConf import*

# Get raw data
HRT   = T1["HRT"]
S1    = T1["S1"]  # [g/L] - COD
XT    = T1["XT"]
S2    = T1["S2"]  # [mmol/L] - VFA
X_1   = T1["X1"]
X_2   = T1["X2"]
Z     = T1["Z"]
C     = T1["C"]
CO2   = T1["CO2"]
B     = T1["B"]
pH    = T1["pH"]
q_C   = T1["q_C"]  # [mmol/L/d]
P_C   = T1["P_C"]  # [atm]
q_M   = T1["q_CH4"]

D     = 1/HRT

S1_in = float(T2["S1in"])    # [gCOD/L]
S2_in = float(T2["S2in"])    # [mmol/L]
C_in  = float(T2["Cin"])    # [mmol/L]
XT_in = float(T2["XTin"])    # [gCOD/L]

# KINETICS

Xr1 = (D*S1, D)
Yr1 = S1

def func1(X,a,b):
    x1,x2 = X
    return a*x1 +a*b*x2 + b*C_d[0]/(1-C_d[0])

#def fun1(beta,X,Y):
    #x1,x2 = X
    #return beta[0]*x1 + beta[0]*beta[1]*x2 + beta[1]*C_d[0]/(1-C_d[0])


popt1, pcov1 = curve_fit(func1, Xr1, Yr1)
#x0 = [0.3,3.5]
#res_lsq1     = least_squares(fun1, x0, args=(Xr1, Yr1))
#print(f"popt: {popt1} lsq: {res_lsq1.}",)
KS1 = popt1[1]                          # [1/d]      Max biomass growth rate
mu1_max = alfa/(1-C_d[0])/popt1[0]      # [g/L]      Half saturation constant

print("KS1:",KS1)
print("mu1_max",mu1_max)

plt.figure(1)
plt.plot(func1(Xr1,popt1[0],popt1[1]),'o',label = 'Estimated')
plt.plot(Yr1,'*', label = 'Real data')
plt.legend()

Xr2 = (D, S2)
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
print(popt2)
print("KS2: %f ; mu2_max: %f; KI2: %f " % ( KS2,  mu2_max, KI2))

plt.figure(2)
plt.plot(func2(Xr2,popt2[0],popt2[1],popt2[2]), 'o',label = 'Estimated')
plt.plot(Yr2,'*', label = 'Real data')
plt.legend()

mu_max = (mu1_max, mu2_max)         # [1/d]      Max biomass growth rate
Ks     = (KS1, KS2)                 # [g/L]      Half saturation constant
KI2    = np.abs(KI2)                # [mmol/L]   Inhibition constant for S2
kd     = np.multiply(C_d,mu_max)


# Hydrolysis

X_hydr = XT
Y_hydr = D*(XT_in-XT)
res_hydr = linregress(X_hydr, Y_hydr)
k_hyd = res_hydr.slope

plt.figure(3)
plt.plot(X_hydr, Y_hydr, 'o', label='original data')
plt.plot(X_hydr, res_hydr.intercept + res_hydr.slope*X_hydr, 'r', label='fitted line')
plt.legend()

# Physical - L/G TRANSFER
pKb = -np.log10(Kb)
fun = 1/(1+10**(pH-pKb))
X3r = C*fun - KH*P_C
Y3r = q_C

mdl3 = linregress(X3r, Y3r)
kLa = mdl3.slope                      # [1/d]      L/G transfer rate
print("kLa:", kLa)
plt.figure(4)
plt.plot(X3r, Y3r, 'o', label = 'data')
plt.plot(X3r, mdl3.intercept + mdl3.slope*X3r, 'r')
plt.legend()

# Yield coefficients: [CODdeg, VFAprof, VFAcons, CO2prod(1), CO2prod(2), CH4prod, hydrolysis]
mu_1 = alfa*D + C_d[0]*mu1_max
mu_2 = alfa*D + C_d[1]*mu2_max


X4r = mu_1*X_1
Y4r = D*(S1_in - S1) + k_hyd*XT
mdl4 = linregress(X4r, Y4r)
k1   = mdl4.slope

X5r = mu_2
Y5r = q_M/X_2
mdl5 = linregress(X5r, Y5r)
X51 = np.array(X5r).reshape((-1,1))


mdl51 = LinearRegression().fit(X51,np.array(Y5r))
k6   = mdl5.slope
k61  = mdl51.coef_

X61 = np.array(D*(S2_in-S2)).reshape((-1,1))
X62 = np.array(D*(S1_in-S1) + k_hyd*XT).reshape((-1,1))
X6r = np.hstack((X61,X62))
Y6r = np.array(q_M)

mdl6 = LinearRegression().fit(X6r,Y6r)
AA = mdl6.coef_[0]
BB = mdl6.coef_[1]/AA

X71 = np.array(D*(S1_in-S1) + k_hyd*XT).reshape((-1,1))
X72 = np.array(q_M).reshape((-1,1))
X7r = np.hstack((X71,X72))
Y7r = np.array(q_C - D*(C_in - C))

# xirls = pylops.optimization.sparsity.IRLS(X7r,Y7r,1)

# print(xirls)


CC1 = 0.2711
DD1 = 0.7934

k2 = BB*k1
k3 = k6/AA
k4 = CC1*k1
k5 = DD1*k6
