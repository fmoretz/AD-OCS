# Modified Version of AMOCO_HN: following the dynamic behaviour from infuent approaching SS
# First try implementation
#

import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from PhysConstants import*
from Influent import*
from Identification import*
from SS_Algebraic import*

# Reactor Configuration

# Substrate definition

# IMPORT FROM ADM1


# PARAMETERS IDENTIFICATION

# SYSTEM

y0 = [y_in[4], 1, 6.5, y_in[2], y_in[0], y_in[1],  y_in[3]]

def f_Model_Dynamics(x,t,alfa,mu_max,Ks,KI2,KH,Pt,kLa,D,y_in,k,kd,N_bac,N_S1):

    XT, X1, X2, Z, S1, S2, C = x


    S1in = y_in[0]
    S2in = y_in[1]
    Cin  = y_in[2]
    Zin  = y_in[3]
    XTin = y_in[4]

    mu1 = mu_max[0]*(S1/(S1+Ks[0]))                     # Monod
    mu2 = mu_max[1]*(S2/(S2+Ks[1]+S2**2/KI2))           # Haldane

    qM  = k[5]*mu2*X2                                   # Methane molar flow [mmol/L/day]
    CO2 = C + S2 - Z                                    # [mmol/L]
    phi = CO2 + KH*Pt + qM/kLa
    Pc  = (phi - (phi**2- 4*KH*Pt*CO2)**0.5)/(2*KH)     # [atm] - Partial pressure CO2
    qC  = kLa*(CO2 - KH*Pc)                             # [mmol/L/d] - Carbon Molar Flow
    CONTROL = KH*Pc**2-Pc*phi + Pt*CO2                  # Control equation, has to be zero

    dXT = D*(XTin - XT) - k[6]*XT                       # Evolution of particulate
    dX1 = (mu1 - alfa*D - kd[0])*X1                     # Evolution of biomass 1 (acidogen.)
    dX2 = (mu2 - alfa*D - kd[1])*X2                     # Evolution of biomass 2 (methanogen)
    dZ  = D*(Zin - Z) + (k[0]*N_S1 - N_bac)*mu1*X1 - N_bac*mu2*X2 + kd[0]*N_bac*X1 + kd[1]*N_bac*X2 # Evolution of alcalinity;
    dS1 = D*(S1in - S1) - k[0]*mu1*X1 + k[6]*XT                                                     # Evolution of organic substrate
    dS2 = D*(S2in - S2) + k[1]*mu1*X1 - k[2]*mu2*X2       # Evolution of VFA
    dC  = D*(Cin - C)   + k[3]*mu1*X1 + k[4]*mu2*X2 - qC  # Evolution of inorganic carbon

    dxdt = [dXT, dX1, dX2, dZ, dS1, dS2, dC]

    return dxdt

tspan = np.linspace(1,50,10000)
YOUT = odeint(f_Model_Dynamics,y0,tspan,args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, y_in, k, kd, N_bac, N_S1))

XT = YOUT[:,0]   # [g/L]
X1 = YOUT[:,1]   # [g/L]
X2 = YOUT[:,2]   # [g/L]
Z  = YOUT[:,3]   # [mmol/L]
S1 = YOUT[:,4]   # [g/L]
S2 = YOUT[:,5]   # [mmol/L]
C  = YOUT[:,6]   # [mmol/L]

# Solver Output
mu1 = np.empty(len(XT))
mu2 = np.empty(len(XT))
CO2 = np.empty(len(XT))
B   = np.empty(len(XT))
phi = np.empty(len(XT))
p_C = np.empty(len(XT))
q_C = np.empty(len(XT))
q_M = np.empty(len(XT))
pH  = np.empty(len(XT))

for x in range(len(XT)):
    mu1[x] = mu_max[0]*(S1[x]/(S1[x]+Ks[0]))                     # Monod
    mu2[x] = mu_max[1]*(S2[x]/(S2[x]+Ks[1]+S2[x]**2/KI2))        # Haldane
    CO2[x] = C[x] + S2[x] - Z[x]
    B[x]   = Z[x] - S2[x]
    phi[x] = CO2[x] + KH*Pt + k[5]/kLa*mu2[x]*X2[x]
    p_C[x]  = (phi[x] - (phi[x]**2- 4*KH*Pt*CO2[x])**0.5)/(2*KH)
    q_C[x] = kLa*(CO2[x] - KH*p_C[x])
    q_M[x] = k[5]*mu2[x]*X2[x]
    pH[x]  = np.real(-np.log10(Kb*CO2[x]/B[x]))

plt.figure(1)
ax1 = plt.subplot()
ax1.plot(tspan, X1, label="Acidogenics")
ax1.plot(tspan, X2, label="Methanogens")
ax1.set_xlabel('time [d]')
ax1.set_ylabel('Microbial Concentration [g/L]')
ax1.grid(True)
ax1.legend()

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('pH Value [-]', color = color)
ax2.plot(tspan, pH, linestyle='dashed',color =color, label="pH")
ax2.tick_params(axis='y', labelcolor=color)

plt.tight_layout()

plt.figure(2)
sub1 = plt.subplot(2,2,1)
sub1.plot(tspan,q_M,label="CH4")
sub1.plot(tspan,q_C,label="CO2")
sub1.set_ylabel('Gas Flowrate [mmol/L/d]')
sub1.set_xlabel('Time [d]')
sub1.grid(True)
sub1.legend()

sub2 = plt.subplot(2,2,2)
sub2.plot(tspan,Z,label="Alkalinity")
sub2.plot(tspan,C,label="In. Carbon")
sub2.set_ylabel('Inorganics Conc. [mmol/L]')
sub2.set_xlabel('Time [d]')
sub2.grid(True)
sub2.legend()

sub3 = plt.subplot(2,2,3)
sub3.plot(tspan,XT,label="Particulate (XT)")
sub3.plot(tspan,S1,label="COD (S1)")
sub3.plot(tspan,S2,label="VFA (S2)")
sub3.set_ylabel('Substrates Conc. [g/L]')
sub3.set_xlabel('Time [d]')
sub3.grid(True)
sub3.legend()

sub4 = plt.subplot(2,2,4)
sub4.plot(tspan,S1,label="COD (S1)")
sub4.plot(tspan,S2,label="VFA (S2) [mmol/L]")
sub4.set_ylabel('Substrates Conc. ')
sub4.set_xlabel('Time [d]')
sub4.grid(True)
sub4.legend()

plt.show()
