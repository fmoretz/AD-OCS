import numpy as np
from matplotlib import pyplot as plt
import time
from pprint import pprint as pp

d_start = 0         # [h] - Start time
d_end   = 5       # [h] - End time
hours   = .01      # [h] - Discretization time

n_times = int((d_end-d_start)/hours)+1  # Number of time steps

print('***Intervals of {hours} hours. \n {n_times} time steps***'.format(hours=hours, n_times=n_times))

t = np.linspace(d_start,d_end,n_times) # time span
R_C = 1.8 # gS/gO2 reaction coefficient
P = 1 # atm
T = 35+273.15 # K
R = 0.0821 # L * atm/mol/K
Vgas = 400#m3
RC = 1.8                                                        # [gS/gO2]     - Stoichiometry of the reaction
alfa = 1                                                        # [-] - Power for sulfur concentration
beta = 1                                                     # [-] - Power for oxygen concentration
k= 0.6                                                     # [1/m3/h*(gS*gO)^proper exp] - Reaction rate constant
tau_head = 4.3                                                  # [h] - Head space residence time
n_in = {
    'H2S':4.54*np.ones(len(t)),
    'O2': 2.54*np.ones(len(t)),
    'SX':  0.000,
    'H2O': 178
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

n_out['H2S'][0] = n_in['H2S'][0]
n_out['O2'][0]  = n_in['O2'][0]
n_out['H2O'][0] = n_in['H2O']

x_out['H2S'][0] = n_in['H2S'][0] / (n_in['H2S'][0] + n_in['O2'][0] + n_in['H2O'] + 2537.39 + 706.49)
x_out['O2'][0]  = n_in['O2'][0]  / (n_in['H2S'][0] + n_in['O2'][0] + n_in['H2O'] + 2537.39 + 706.49)

c_out['H2S'][0] = x_out['H2S'][0] * P/R/T # mol/L
c_out['O2'][0]  = x_out['O2'][0]  * P/R/T # mol/L

m_out['H2S'][0] = c_out['H2S'][0] * PM['H2S'] * 1000 # g/m3
m_out['O2'][0]  = c_out['O2'][0]  * PM['O2'] * 1000 # g/m3

print('Simulation is on going...')
t0 = time.time()
r_SOB[0] = max(0,(k * m_out['H2S'][0]**alfa * m_out['O2'][0]**beta)) # gS/m3/h
for i in range(1, len(t)):
#                                                                                     CH4       CO2
    x_out['H2S'][i-1] = n_out['H2S'][i-1] / (n_out['H2O'][i-1] + n_out['H2S'][i-1] + n_out['O2'][i-1] + 2537.39 + 706.49)
    x_out['O2'][i-1]  = n_out['O2'][i-1]  / (n_out['H2O'][i-1] + n_out['H2S'][i-1] + n_out['O2'][i-1] + 2537.39 + 706.49)
    
    c_out['H2S'][i-1] = x_out['H2S'][i-1] * P/R/T # mol/L
    c_out['O2'][i-1]  = x_out['O2'][i-1]  * P/R/T # mol/L
    
    m_out['H2S'][i-1] = c_out['H2S'][i-1] * PM['H2S'] * 1000 # g/m3
    m_out['O2'][i-1]  = c_out['O2'][i-1]  * PM['O2'] * 1000 # g/m3
    

    r_SOB[i] = (k * m_out['H2S'][i-1]**alfa * m_out['O2'][i-1]**beta)

    n_out['H2S'][i] = max(1e-16,(n_in['H2S'][i]*dt + (n_out['H2S'][i-1])  - r_SOB[i]*dt* Vgas / PM['H2S'])/(1+dt/tau_head))
    n_out['O2'][i]  = max(1e-16,(n_in['O2'][i]*dt  + (n_out['O2'][i-1])  - r_SOB[i]*dt* Vgas / (R_C * PM['O2']))/(1+dt/tau_head))
    n_out['SX'][i]  = (n_in['SX']*dt  + (n_out['SX'][i-1])  + r_SOB[i]* dt* Vgas / PM['SX'])/(1+dt/tau_head)
    
    n_out['H2O'][i] = n_in['H2O'] + (n_out['H2S'][i-1] - n_in['H2S'][i]) # n_in['H2O'] + nu['H2O'] * (k * m_out['H2S'][i]**alfa * m_out['O2'][i]**beta)*( dt/(1+dt) ) * Vgas /



Efficiency = (1-(n_out['H2S'][:]/tau_head/n_in['H2S']))*100
    
print(n_out['O2'][:])
print('Simulation ended')
print('time for procedure: ', time.time() - t0, 'sec')

plt.figure()

plt.plot(t, n_out['H2S'], label = 'H2S', linewidth=2)
plt.plot(t, n_out['O2'], '--', label = 'O2',  linewidth=2)
plt.plot(t, n_out['SX'], '-.', label = 'SX',  linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('mol')
plt.legend(loc='best')


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


