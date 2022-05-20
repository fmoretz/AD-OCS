# import numerical libraries
import numpy as np
from AMOCOHN_main import*
import matplotlib.pyplot as plt
from PhysConstants import*

# define logistic function
def logi(x,a,b,c):
    return a/(1+np.exp(-b*(x-c)))

# define gompertz function
def gompertz(x,a,b,c):
    return a*np.exp(-np.exp(b*np.exp(1)/a*(c-x) - 1))

# define lokta-volterra function
def lokta_volterra(x,a,b,c,d):
    return a*np.exp(b*x) + c*np.exp(d*x)

# retrieve data from amocohn_main
for x in range(len(XT)):
    mu1[x]  = mu_max[0]*(S1[x]/(S1[x]+Ks[0]))                     # Monod
    mu2[x]  = mu_max[1]*(S2[x]/(S2[x]+Ks[1]+S2[x]**2/KI2))        # Haldane
    CO2[x]  = C[x] + S2[x] - Z[x]
    B[x]    = Z[x] - S2[x]
    phi[x]  = CO2[x] + KH*Pt + k[5]/kLa*mu2[x]*X2[x]
    p_C[x]  = (phi[x] - (phi[x]**2- 4*KH*Pt*CO2[x])**0.5)/(2*KH)
    q_C[x]  = kLa*(CO2[x] - KH*p_C[x])
    q_M[x]  = k[5]*mu2[x]*X2[x]
    pH[x]   = np.real(-np.log10(Kb*CO2[x]/B[x]))
    
q_tot = q_C+q_M
x_M   = np.divide(q_M,q_tot)

# Sulfur production model: dSs/dt = (1-y)/y * dXs/dt
Ss_max = 0.03*XT_in[0]    # maximum sulfur concentration ====
Lam    = mu2_max         # Sulfur production lag-time, eqaul to maximum S2 consumption rate [d] ===
Xs_max = 0.035*S2[-1]*k[2]/MW_AAc    # maximum sulfur concentration g/L ====

print('\nSulfur production model: dSs/dt = (1-y)/y * dXs/dt')
print(f'Max sulfur concentration: {Ss_max}')
print(f'Max SRB consumption rate: {Xs_max}')
print(f'Sulfur production lag-time: {Lam}')

Ss_logi = np.empty(len(XT))
Xs_logi = np.empty(len(XT))
Ss_gomp = np.empty(len(XT))
Xs_gomp = np.empty(len(XT))

DS1 = np.empty(len(XT))
DS2 = np.empty(len(XT))
DXT = np.empty(len(XT))
DX1 = np.empty(len(XT))
DX2 = np.empty(len(XT))

for i in range(len(XT)):
    # Species differences
    DS1[i] = - S1[0] + S1[i]
    DS2[i] = - S2[0] + S2[i]
    DX1[i] = - X1[0] + X1[i]
    DXT[i] = - XT[0] + XT[i]
    DX2[i] = - X2[0] + X2[i]
    
    # Logistic function
    Ss_logi[i] = logi(tspan[i], Ss_max, DX2[i], Lam)
    Xs_logi[i] = logi(tspan[i], Xs_max, DX2[i], Lam)

    # Gompertz function
    Ss_gomp[i] = gompertz(tspan[i], Ss_max, DX2[i], Lam)
    Xs_gomp[i] = gompertz(tspan[i], Xs_max, DX2[i], Lam)
    

plt.figure(100)
plt.subplot(2,1,1)
plt.plot(tspan,Ss_logi, label="Logistic")
plt.plot(tspan,Ss_gomp, label="Gompertz")
plt.ylabel('Sulfur Concentration g/L')
plt.grid(True); plt.legend(loc='best')

plt.subplot(2,1,2)
plt.plot(tspan,Xs_logi, label="Logistic")
plt.plot(tspan,Xs_gomp, label="Gompertz")
plt.ylabel('Biomass Concentration g/L')
plt.xlabel('Time [d]')
plt.grid(True); plt.legend(loc='best')

plt.figure(101)
plt.subplot(5,1,1)
plt.plot(tspan, DS1, label="S1")
plt.grid(True)
plt.legend()

plt.subplot(5,1,2)
plt.plot(tspan, DS2, label="S2")
plt.grid(True)
plt.legend()

plt.subplot(5,1,3)
plt.plot(tspan, DXT, label="XT")
plt.ylabel('Difference Concentration [g/L]')
plt.grid(True)
plt.legend()

plt.subplot(5,1,4)
plt.plot(tspan, DX2, label="X2")
plt.grid(True)
plt.legend()

plt.subplot(5,1,5)
plt.plot(tspan, DX1, label="X1")
plt.xlabel('Time [d]')
plt.grid(True)
plt.legend()

plt.show()