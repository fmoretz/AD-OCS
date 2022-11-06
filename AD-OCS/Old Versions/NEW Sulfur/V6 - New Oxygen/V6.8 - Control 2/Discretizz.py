import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import time
from pprint import pprint as pp

def SS_evaluation(F_out, Fin, tau, Vgas):
    global alpha, beta
    F_out_H2S = F_out[0]
    F_out_O2 = F_out[1]
    F_out_H2O = F_out[2]

    x_out = {'CH4': Fin['CH4']/(Fin['CH4']+Fin['CO2']+sum(F_out)),
            'CO2': Fin['CO2']/(Fin['CH4']+Fin['CO2']+sum(F_out)),
            'H2S': F_out_H2S/(Fin['CH4']+Fin['CO2']+sum(F_out)),
            'H2O': F_out_H2O/(Fin['CH4']+Fin['CO2']+sum(F_out)),
            'O2': F_out_O2/(Fin['CH4']+Fin['CO2']+sum(F_out)),
    }
    c_out = (x_out[comp][i]*P/R/T for comp in x_out.keys())
    w_out = (c_out[comp]*PM[comp] for comp in c_out.keys() if comp != 'CH4' and comp != 'CO2')

    r_SOB = k*w_out['H2S']**alpha+w_out['H2O']**beta

    eq1 = (F_in['H2S']-F_out_H2S + r_SOB/PM['H2S']*Vgas)
    eq2 = (F_in['O2']-F_out_O2 + r_SOB*Vgas)
    eq3 = (F_in['H2O']-F_out_H2O + r_SOB*Vgas)
    return [eq1, eq2, eq3]

d_start = 0         # [h] - Start time
d_end   = 60       # [h] - End time
hours   = .01      # [h] - Discretization time

n_times = int((d_end-d_start)/hours)+1  # Number of time steps

print('***Intervals of {hours} hours. \n {n_times} time steps***'.format(hours=hours, n_times=n_times))

t = np.linspace(d_start,d_end,n_times) # time span
R_C = 1.8 # gS/gO2 reaction coefficient
P = 1 # atm
T = 35+273.15 # K
P_norm = 1 # atm
T_norm = 20+273.15 # K
R = 0.0821 # L * atm/mol/K
Rgas_m3_atm_K = 0.08205746e-3 # m^3*atm/(mol*K)
Vgas = 350*np.ones(len(t))#m3
RC = 1.8                                                        # [gS/gO2]     - Stoichiometry of the reaction
alfa = 1                                                        # [-] - Power for sulfur concentration
beta = 0.2                                                      # [-] - Power for oxygen concentration
k= 0.6                                                          # [1/m3/h*(gS*gO)^proper exp] - Reaction rate constant
tau_head = 4.3*np.ones(len(t))                                                  # [h] - Head space residence time
F_in = {
    'CH4': 2537.39*np.ones(len(t)),
    'CO2':  706.49*np.ones(len(t)),
    'H2S': 4.54*np.ones(len(t)),
    'H2O': 178*np.ones(len(t)),
    'O2':  2.2*np.ones(len(t)),
    'SX':  0.000,    
    'total': np.zeros(len(t))   
} # mol/hours

PM = {
    'H2S': 34.1,
    'O2':  36.0,
    'SX':  32.1,
    'H2O': 18.0
}

dt = (t[-1] - t[0])/len(t)

n_out = {
    'CH4': np.zeros(len(t)),
    'CO2': np.zeros(len(t)),
    'H2S': np.zeros(len(t)),
    'H2O': np.zeros(len(t)),
    'O2':  np.zeros(len(t)),
    'SX':  np.zeros(len(t)),    
} # mol
F_out = {
    'CH4': np.zeros(len(t)),
    'CO2': np.zeros(len(t)),
    'H2S': np.zeros(len(t)),
    'H2O': np.zeros(len(t)),
    'O2':  np.zeros(len(t)),
    'SX':  np.zeros(len(t)),
    'total':np.zeros(len(t))   
} # mol/hours
Q_out = {
    'CH4': np.zeros(len(t)),
    'CO2': np.zeros(len(t)),
    'H2S': np.zeros(len(t)),
    'H2O': np.zeros(len(t)),
    'O2':  np.zeros(len(t)),
    'SX':  np.zeros(len(t)),
    'total':np.zeros(len(t))   
} # mol/hours

c_out = {
    'H2S': np.zeros(len(t)),
    'O2':  np.zeros(len(t))
} # mol/m3

w_out = {
    'H2S': np.zeros(len(t)),
    'O2':  np.zeros(len(t))
} # g/m3

x_out = {
    'H2S': np.zeros(len(t)),
    'O2':  np.zeros(len(t))
} # [-]

r_SOB = np.zeros(len(t))

F_out['H2S'][0] = 1.5865
F_out['O2'][0]  = 0.5
F_out['H2O'][0] = 180.9
F_out['CH4'][0] = F_in['CH4'][0]
F_out['CO2'][0] = F_in['CO2'][0]


n_out['H2S'][0] = F_out['H2S'][0]*tau_head[0]
n_out['O2'][0]  = F_out['O2'][0]*tau_head[0]
n_out['H2O'][0] = F_out['H2O'][0]*tau_head[0]
n_out['SX'][0] =  13.5

digester_out = pd.DataFrame()
print('Simulation is on going...')
t0 = time.time()
# r_SOB[0] = (k * w_out['H2S'][0]**alfa * w_out['O2'][0]**beta)) # gS/m3/h
H2S_treshold = 400 # ppm

for i in range(1, len(t)):#   
    print('Time step: {t} h'.format(t=t[i], i=i))
    it = 1
    while True:         
            
        F_in['total'][i] = sum(F_in[comp][i] for comp in F_in.keys() if comp != 'SX')
        if i > len(t)/3: # and i < len(t)/2:
            F_in['H2S'][i] = 1.3*F_in['H2S'][0]   
            
        x_out['H2S'][i-1] = F_out['H2S'][i-1]*tau_head[i-1] / (F_out['H2O'][i-1]*tau_head[i-1] + F_out['H2S'][i-1]*tau_head[i-1] + F_out['O2'][i-1] *tau_head[i-1]+  F_out['CH4'][i-1]*tau_head[i-1] +  F_out['CO2'][i-1]*tau_head[i-1])
        x_out['O2'] [i-1]  = F_out['O2'][i-1]*tau_head[i-1] / (F_out['H2O'][i-1]*tau_head[i-1] + F_out['H2S'][i-1]*tau_head[i-1] + F_out['O2'][i-1] *tau_head[i-1]+  F_out['CH4'][i-1]*tau_head[i-1] +  F_out['CO2'][i-1]*tau_head[i-1])
        
        c_out['H2S'][i-1] = x_out['H2S'][i-1] * P/R/T # mol/L
        c_out['O2'][i-1]  = x_out['O2'][i-1]  * P/R/T # mol/L
        
        w_out['H2S'][i-1] = c_out['H2S'][i-1] * PM['H2S'] * 1000 # g/m3
        w_out['O2'][i-1]  = c_out['O2'][i-1]  * PM['O2'] * 1000 # g/m3
        
        r_SOB[i] = (k * w_out['H2S'][i-1]**alfa * w_out['O2'][i-1]**beta) # g/m3/h

        n_out['H2S'][i] = max(1e-16,(F_in['H2S'][i]*dt + (n_out['H2S'][i-1]) - (r_SOB[i]*dt* Vgas[i] / PM['H2S']))/(1+dt/tau_head[i]))
        n_out['O2'][i]  = max(1e-16,(F_in['O2'][i]*dt  + (n_out['O2'][i-1])  - (r_SOB[i]*dt* Vgas[i] / (R_C * PM['O2'])))/(1+dt/tau_head[i]))
        n_out['SX'][i]  = (F_in['SX']*dt  + (n_out['SX'][i-1])  + (r_SOB[i]* dt* Vgas[i]) / PM['SX'])/(1+dt/tau_head[i])
        
        n_out['H2O'][i] = (F_in['H2O'][i]*dt + (n_out['H2O'][i-1])  + (r_SOB[i]*dt* Vgas[i] / PM['H2S']))/(1+dt/tau_head[i]) # F_in['H2O'] + nu['H2O'] * (k * w_out['H2S'][i]**alfa * w_out['O2'][i]**beta)*( dt/(1+dt) ) * Vgas /

        F_out['H2S'][i] = n_out['H2S'][i]/tau_head[i]
        F_out['O2'][i]  = n_out['O2'][i]/tau_head[i]
        
        F_out['H2O'][i] = n_out['H2O'][i]/tau_head[i]
        F_out['CO2'][i] = F_in['CO2'][i]
        F_out['CH4'][i] = F_in['CH4'][i]
        F_out['total'][i] = F_out['H2S'][i] + F_out['O2'][i] + F_out['H2O'][i] + F_out['CO2'][i] + F_out['CH4'][i]
        dig_iter = pd.DataFrame()
        dig_iter = pd.DataFrame({'t': t[i], 
                                't_out': t[i]+tau_head[i], 
                                'F_CH4': F_out['CH4'][i],
                                'F_CO2': F_out['CO2'][i],
                                'F_H2S': F_out['H2S'][i],
                                'F_O2' : F_out['O2'][i],
                                'F_H2O': F_out['H2O'][i],
                                'F_tot': F_out['total'][i],
                                'x_CH4': F_out['CH4'][i]/F_out['total'][i],
                                'x_CO2': F_out['CO2'][i]/F_out['total'][i],
                                'x_H2S': F_out['H2S'][i]/F_out['total'][i],
                                'x_H2O': F_out['H2O'][i]/F_out['total'][i],
                                'x_O2' : F_out['O2'][i]/F_out['total'][i],
                                'Q_CH4': F_out['CH4'][i]*T_norm*Rgas_m3_atm_K/P_norm,
                                'Q_CO2': F_out['CO2'][i]*T_norm*Rgas_m3_atm_K/P_norm,
                                'Q_H2S': F_out['H2S'][i]*T_norm*Rgas_m3_atm_K/P_norm,
                                'Q_H2O': F_out['H2O'][i]*T_norm*Rgas_m3_atm_K/P_norm,
                                'Q_O2' : F_out['O2'][i]*T_norm*Rgas_m3_atm_K/P_norm,
                                'Q_tot': F_out['total'][i]*T_norm*Rgas_m3_atm_K/P_norm,
                                }, index=[i])

        if dig_iter['x_H2S'][i] < H2S_treshold:
            break
        else:
            print('H2S above treshold, time step: ', i)
            F_in['O2'][i-1] = F_in['O2'][i-1]*1.05

    digester_out = pd.concat([digester_out, dig_iter])



                            
                          
print(digester_out)                          
    
Efficiency = (1-(n_out['H2S'][:]/tau_head/F_in['H2S']))*100
    
print('Last F_out H2S = ', F_out['H2S'][-1], 'mol/h')
print('Last F_out O2  = ', F_out['O2'][-1], 'mol/h')
print('Last F_out H2O = ', F_out['H2O'][-1], 'mol/h')
print('Last F_out CO2 = ', F_out['CO2'][-1], 'mol/h')
print('Last n_out SX  = ', n_out['SX'][-1], 'mol')
print('Simulation ended')
print('time for procedure: ', time.time() - t0, 'sec')

plt.figure()
plt.subplots_adjust(hspace=0.7)
plt.subplot(3,1,1)
'CSTR composition'
plt.plot(t, n_out['H2S'], label = 'H2S', linewidth=2)
plt.plot(t, n_out['O2'], '--', label = 'O2',  linewidth=2)
plt.plot(t, n_out['SX'], '-.', label = 'SX',  linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('mol')
plt.legend(loc='best')
plt.subplot(3,1,2)
plt.title('H2S Flowrate')
plt.plot(t, F_out['H2S'], label = 'out', linewidth=2)
plt.plot(t, F_in['H2S'], '--', label = 'in',  linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('mol/h'); plt.legend(loc='best')
plt.subplot(3,1,3)
plt.title('O2Flowrate')
plt.plot(t, F_out['O2'], label = 'out', linewidth=2)
plt.plot(t, F_in['O2'], '--', label = 'in',  linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('mol/h'); plt.legend(loc='best')
plt.figure()
plt.title('H2O')
plt.plot(t, F_out['H2O'], label = 'out', linewidth=2)
plt.plot(t, F_in['H2O'], '--', label = 'in',linewidth=2)
plt.legend()
plt.figure()
plt.subplot(2,1,1)
plt.title('Reaction Rate')
plt.plot(t[1:], r_SOB[1:]*24, label = 'r_SOB', linewidth=2)
plt.plot(t, np.array(np.mean(r_SOB*24))*np.ones(len(t)), '--', label = 'mean',  linewidth=2)
plt.legend(loc='best')
plt.grid(True); plt.xlabel('hours'); plt.ylabel('mgS/L/d')
plt.subplot(2,1,2)
plt.plot(t, Efficiency, label = 'Efficiency', linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('Efficiency, %')
print('mean r_SOB: ', np.mean(r_SOB*24), 'mgS/L/d')

plt.figure()
plt.title('Digester Out')
plt.plot(digester_out['t'], digester_out['x_H2S']*1e+6, label = 'H2S', linewidth=2)
plt.ylabel(ylabel='x_H2S, ppm')
plt.show()


