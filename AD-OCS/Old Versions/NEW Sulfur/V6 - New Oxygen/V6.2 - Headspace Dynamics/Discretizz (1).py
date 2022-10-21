import numpy as np
from matplotlib import pyplot as plt
import time
from pprint import pprint as pp

k    = 0.6
alfa = 1
beta = 0.1
R_C = 2 # gS/gO2 reaction coefficient
P = 1 # atm
T = 35+273.15 # K
R = 0.0821 # L * atm/mol/K
Vgas = 350 #m3

t = np.linspace(0, 7, 100000)    # hours

n_in = {
    'H2S': 5.153,
    'O2':  5.153/2,
    'SX':  0.000,
    'H2O': 191.3
} # mol/hours

PM = {
    'H2S': 34.1,
    'O2':  36.0,
    'SX':  32.1,
    'H2O': 18.0
}


dt = (t[-1] - t[0])/len(t)

n_out = {
    'H2S': np.zeros(len(t)),
    'O2':  np.zeros(len(t)),
    'SX':  np.zeros(len(t)),
    'H2O': np.zeros(len(t))
} # mol/hours

c_out = {
    'H2S': np.zeros(len(t)),
    'O2':  np.zeros(len(t))
} # mol/hours

m_out = {
    'H2S': np.zeros(len(t)),
    'O2':  np.zeros(len(t))
} # mol/hours

x_out = {
    'H2S': np.zeros(len(t)),
    'O2':  np.zeros(len(t))
} # mol/hours

r_SOB = np.zeros(len(t))

n_out['H2S'][0] = n_in['H2S']
n_out['O2'][0]  = n_in['O2']
n_out['H2O'][0] = n_in['H2O']

x_out['H2S'][0] = n_in['H2S'] / (sum(n_in.values()) - n_in['SX'] + 2537.39 + 706.49)
x_out['O2'][0]  = n_in['O2']  / (sum(n_in.values()) - n_in['SX'] + 2537.39 + 706.49)

c_out['H2S'][0] = x_out['H2S'][0] * P/R/T # mol/L
c_out['O2'][0]  = x_out['O2'][0]  * P/R/T # mol/L

m_out['H2S'][0] = c_out['H2S'][0] * PM['H2S'] * 1000 # g/m3
m_out['O2'][0]  = c_out['O2'][0]  * PM['O2'] * 1000 # g/m3

print('Simulation is on going...')
t0 = time.time()
r_SOB[0] = (k * m_out['H2S'][0]**alfa * m_out['O2'][0]**beta) # gS/m3/h
for i in range(1, len(t)):
#                                                                                     CH4       CO2
    x_out['H2S'][i-1] = n_out['H2S'][i-1] / (n_out['H2O'][i-1] + n_out['H2S'][i-1] + n_out['O2'][i-1] + 2537.39 + 706.49)
    x_out['O2'][i-1]  = n_out['O2'][i-1]  / (n_out['H2O'][i-1] + n_out['H2S'][i-1] + n_out['O2'][i-1] + 2537.39 + 706.49)
    
    c_out['H2S'][i-1] = x_out['H2S'][i-1] * P/R/T # mol/L
    c_out['O2'][i-1]  = x_out['O2'][i-1]  * P/R/T # mol/L
    
    m_out['H2S'][i-1] = c_out['H2S'][i-1] * PM['H2S'] * 1000 # g/m3
    m_out['O2'][i-1]  = c_out['O2'][i-1]  * PM['O2'] * 1000 # g/m3
    

    r_SOB[i] = (k * m_out['H2S'][i-1]**alfa * m_out['O2'][i-1]**beta)

    n_out['H2S'][i] = (n_in['H2S']*dt + (n_out['H2S'][i-1])  - r_SOB[i]*dt* Vgas / PM['H2S'])/(1+dt)
    n_out['O2'][i]  = (n_in['O2']*dt  + (n_out['O2'][i-1])  - r_SOB[i]*dt* Vgas / (R_C * PM['O2']))/(1+dt)
    n_out['SX'][i]  = (n_in['SX']*dt  + (n_out['SX'][i-1])  + r_SOB[i]* dt* Vgas / PM['SX'])/(1+dt)
    
    n_out['H2O'][i] = n_in['H2O'] + (n_out['H2S'][i-1] - n_in['H2S']) # n_in['H2O'] + nu['H2O'] * (k * m_out['H2S'][i]**alfa * m_out['O2'][i]**beta)*( dt/(1+dt) ) * Vgas /


Efficiency = (1-(n_out['H2S'][:]/n_in['H2S']))*100
    

print('Simulation ended')
print('time for procedure: ', time.time() - t0, 'sec')

plt.figure()
plt.subplot(2,1,1)
plt.plot(t, n_out['H2S'], label = 'H2S', linewidth=2)
plt.plot(t, n_out['O2'], '--', label = 'O2',  linewidth=2)
plt.plot(t, n_out['SX'], '-.', label = 'SX',  linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('mol/hours')
plt.legend(loc='best')
plt.subplot(2,1,2)
plt.plot(t, m_out['H2S'], label = 'H2S', linewidth=2)
plt.plot(t, m_out['O2'], '--', label = 'O2',  linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('g/m3')

plt.figure()
plt.subplot(2,1,1)
plt.title('Reaction Rate')
plt.plot(t, r_SOB*24, label = 'r_SOB', linewidth=2)
plt.plot(t, np.array(np.mean(r_SOB*24))*np.ones(len(t)), '--', label = 'mean',  linewidth=2)
plt.legend(loc='best')
plt.grid(True); plt.xlabel('hours'); plt.ylabel('mgS/L/d')
plt.subplot(2,1,2)
plt.plot(t, Efficiency, label = 'Efficiency', linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('Efficiency, %')
print('mean r_SOB: ', np.mean(r_SOB*24), 'mgS/L/d')
plt.show()


