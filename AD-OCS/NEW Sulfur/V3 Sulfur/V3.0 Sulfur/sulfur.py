# import numerical libraries
import numpy as np
from main import*
import matplotlib.pyplot as plt
from PhysConstants import*
from plot import fitsubplotter

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
Y_srb = 0.035*64/1000

Xs = np.empty(len(XT))
lam = np.empty(len(XT))
Ss  = np.empty(len(XT))
y_S  = np.empty(len(XT))
Ss_max= np.empty(len(XT))
Xs_max = np.empty(len(XT))
mu_max_srb = np.empty(len(XT))
rho_srb = np.empty(len(XT))

H_S  = 0.1367719*T + 2.0180702               # [atm]        - Henry's constant CO2 at [T]°C - Partial Pressure Relation: Linear regression on Perry's data
H_S  = 0.00030308*(T**2) + 0.11276471*T +2.44557423
KH_S = 1/(H_S*100/55342)                    # [mmol/L/atm] - Henry's constant CO2 at [T]°C - Concentrations Relation
H_S  = 1/KH_S 


for i in range(len(XT)):
    # Species differences
    Ss_max[i] = 0.02*y_influent[i,4]*1000/64*S2[i]/y_influent[i,4]
    Xs_max[i] = Y_srb/(1-Y_srb)*Ss_max[i]   # maximum sulfur concentration g/L
    
    
    if i < 1:
        mu_max_srb[i] = 0
        lam[i]    = 0
        Xs [i] =gompertz(t_span[i], Xs_max[i], mu_max_srb[i], lam[i])
        rho_srb[i] = growth_SRB(t_span[i], Xs_max[i], mu_max_srb[i], lam[i])
    else:
            # Gompertz function
        mu_max_srb[i] = max(1e-16,(- X2[0] + X2[i]))/(t_span[i] - t_span[0])
        lam[i]    = 0
        Xs[i]     = gompertz(t_span[i], Xs_max[i], mu_max_srb[i], lam[i])
        rho_srb[i] = growth_SRB(t_span[i], Xs_max[i], mu_max_srb[i], lam[i])

    Ss[i]     = (1-Y_srb)/(Y_srb)*Xs[i]
    y_S[i]    = (H_S*Ss[i])/1

print(mu_max_srb)

plt.figure(100)
plt.subplot(4,1,1)
plt.plot(t_span, Xs, label="Gompertz")
plt.ylabel('Xs [g/L')
plt.grid(True)

plt.subplot(4,1,2)
plt.grid(True)
plt.plot(t_span, rho_srb)

# plt.plot(t_span, mu_max_srb,'--')
plt.xlabel('Time [d]')
plt.ylabel('X growt [g/L/d]')

plt.subplot(4,1,3)
plt.grid(True)
plt.plot(t_span, Ss)
plt.xlabel('Time [d]')
plt.ylabel('Ss [mmol/L]')

plt.subplot(4,1,4)
plt.grid(True)
plt.plot(t_span, y_S)
plt.xlabel('Time [d]')
plt.ylabel('yS [-]')

plt.show()