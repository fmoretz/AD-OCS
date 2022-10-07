# import numerical libraries
import numpy as np
from main_V3 import*
import matplotlib.pyplot as plt
from PhysConstants import*
#from plot import fitsubplotter

# define gompertz function
def gompertz(x,a,b,c):
    return a*np.exp(-np.exp(b*np.exp(1)/a*(c-x) - 1))
def growth_SRB(x,a,b,c):
    return b*np.exp((np.exp(1)*b/a*(c-x))-np.exp(np.exp(1)*b/a*(c-x)+1)+2)


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
Xs_max = 1    # maximum sulfur concentration g/L

Xs = np.empty(len(XT))
lam = np.empty(len(XT))
mu_max_srb = np.empty(len(XT))
mu_srb = np.empty(len(XT))
for i in range(len(XT)):
    # Species differences
    
    mu_max_srb[i] = (- X2[0] + X2[i])/(tspan[i] - tspan[0])
    lam[i] = mu_max_srb[i] - Xs_max/(mu_max_srb[i]*np.exp(1))
    
    print(f'lam[{i}] = {lam[i]} - {mu_max_srb[i]}')

    # Gompertz function
    Xs[i] = gompertz(tspan[i], Xs_max, mu_max_srb[i], lam[i])
    mu_srb[i] = growth_SRB(tspan[i], Xs_max, mu_max_srb[i], lam[i])

    

plt.figure(100)
plt.subplot(3,1,1)
plt.plot(tspan, Xs, label="Gompertz")
plt.ylabel('Biomass Concentration g/L')
plt.grid(True)

plt.subplot(3,1,2)
plt.grid(True)
plt.plot(tspan, mu_srb)
plt.xlabel('Time [d]')
plt.ylabel('X growt [g/L/d]')

plt.subplot(3,1,3)
plt.grid(True)
plt.plot(tspan, mu_srb/mu_max_srb)
plt.xlabel('Time [d]')
plt.ylabel('Lag-Time [d]')
print(len(XT))
plt.show()