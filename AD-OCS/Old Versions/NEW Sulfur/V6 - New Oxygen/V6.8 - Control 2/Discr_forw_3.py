import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import time
from pprint import pprint as pp
'''  *** Development of discr_forw_1 *** ATTENTION: measure units, x_H2S are mol/mol; threshold is ppm. If any change is needed take into accont.
Control algorithm for the oxygen injection: it is based on increasing the ONLY THE PREVIOUS ([i-1]) influent flow rate
in order to reach the desired H2S concentration in the effluent.
- V2: Limit to maximum influent flow rate
- V3: Going back in time'''
P_norm = 1 # atm
T_norm = 20+273.15 # K
R = 0.0821 # L * atm/mol/K
Rgas_m3_atm_K = 0.08205746e-3 # m^3*atm/(mol*K)

t_change = 24
d_start = 0         # [h] - Start time
d_end   = 60       # [h] - End time
hours   = .1      # [h] - Discretization time

n_times = int((d_end-d_start)/hours)+1  # Number of time steps

air = 2 # 1 for air, 2 for O2

H2S_treshold = 450 # ppm
Q_in_max = 10 # m3/h air or O2 according to the "air" variable
N_in_max = 10 # mol/h
Q_in = 0.2                                                 # [m3/h] - Inlet volumetric flow - Gas (air or oxygen) defined in the input section
N_in_O2 = 2.2

control = 2                                             # 1 for control on Qin, 2 for control on N_in_O2

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

if control == 1:        
    if air == 1:
        N_in_O2 = Q_in*0.21*P_norm/(Rgas_m3_atm_K*(T_norm+273.15))    # [mol/h] - Inlet molar flow - Oxygen
        N_in_N2 = Q_in*0.79*P_norm/(Rgas_m3_atm_K*(T_norm+273.15))    # [mol/h] - Inlet molar flow - Nitrogen
    elif air == 2:
        N_in_O2 = Q_in*P_norm/(Rgas_m3_atm_K*(T_norm+273.15))         # [mol/h] - Inlet molar flow - Oxygen
        N_in_N2 = 0.000
    else:
        print('Error: air variable not defined')
    F_in_max = Q_in_max*0.21*P_norm/(Rgas_m3_atm_K*(T_norm+273.15))         # [mol/h] - Inlet molar flow - Oxygen
elif control == 2:
    N_in_O2 = N_in_O2
    N_in_N2 = 0.000
    F_in_max = N_in_max

F_in = {
        'CH4': 2537.39*np.ones(len(t)),
    'CO2':  706.49*np.ones(len(t)),
    'H2S': 4.54*np.ones(len(t)),
    'H2O': 178*np.ones(len(t)),
    'O2':  N_in_O2*np.ones(len(t)),
    'SX':  0.000*np.ones(len(t)), 
    'N2':  N_in_N2*np.ones(len(t)),  
    'total': 0.000*np.ones(len(t))
   } # mol/hours
for i in range(len(t)):
    F_in['total'][i] = sum(F_in[comp][i] for comp in F_in.keys() if comp != 'total')

PM = {
    'H2S': 34.1,
    'O2':  36.0,
    'SX':  32.1,
    'H2O': 18.0
}

dt = (t[-1] - t[0])/len(t)

n_mol = {
    'CH4': np.zeros(len(t)),
    'CO2': np.zeros(len(t)),
    'H2S': np.zeros(len(t)),
    'H2O': np.zeros(len(t)),
    'O2':  np.zeros(len(t)),   
    'SX':  np.zeros(len(t)),
    'N2':  np.zeros(len(t)),   
} # mol
F_out = {
    'CH4': np.zeros(len(t)),
    'CO2': np.zeros(len(t)),
    'H2S': np.zeros(len(t)),
    'H2O': np.zeros(len(t)),
    'O2':  np.zeros(len(t)),
    'SX':  np.zeros(len(t)),
    'N2':  np.zeros(len(t)),
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

F_out['H2S'][0] = 1.51
F_out['O2'][0]  = 0.7
F_out['H2O'][0] = 181.03
F_out['CH4'][0] = F_in['CH4'][0]
F_out['CO2'][0] = F_in['CO2'][0]

n_mol['H2S'][0] = F_out['H2S'][0]*tau_head[0]
n_mol['O2'][0]  = F_out['O2'][0]*tau_head[0]
n_mol['H2O'][0] = F_out['H2O'][0]*tau_head[0]
n_mol['SX'][0] =  13.05
# r_SOB[0] = (k * w_out['H2S'][0]**alfa * w_out['O2'][0]**beta)) # gS/m3/h

print('Simulation is on going...')
t0 = time.time()

for i in range(0, len(t)-1):#   
    if t[i] > t_change: # and i < len(t)/2:
        F_in['H2S'][i] = 1.3*F_in['H2S'][0]   
    
    while True: 

        x_out['H2S'][i]  = F_out['H2S'][i]*tau_head[i] / (F_out['H2O'][i]*tau_head[i]  + F_out['H2S'][i]*tau_head[i] + F_out['O2'][i] *tau_head[i]+  F_out['CH4'][i]*tau_head[i] +  F_out['CO2'][i]*tau_head[i] + F_out['N2'][i]*tau_head[i])
        x_out['O2'] [i]  = F_out['O2'][i]*tau_head[i] /  (F_out['H2O'][i]*tau_head[i]  + F_out['H2S'][i]*tau_head[i] + F_out['O2'][i] *tau_head[i]+  F_out['CH4'][i]*tau_head[i] +  F_out['CO2'][i]*tau_head[i] + F_out['N2'][i]*tau_head[i])
        
        c_out['H2S'][i] = x_out['H2S'][i] * P/R/T # mol/L
        c_out['O2'][i]  = x_out['O2'][i]  * P/R/T # mol/L
        
        w_out['H2S'][i] = c_out['H2S'][i] * PM['H2S'] * 1000 # g/m3
        w_out['O2'][i]  = c_out['O2'][i]  * PM['O2'] * 1000 # g/m3
        
        r_SOB[i] = (k * w_out['H2S'][i]**alfa * w_out['O2'][i]**beta) # gS/m3/h = mgS/L/h

        n_mol['H2S'][i+1] = n_mol['H2S'][i] + (F_in['H2S'][i] - n_mol['H2S'][i]/tau_head[i] - r_SOB[i]*Vgas[i]/PM['SX'])*dt
        n_mol['O2'][i+1]  = n_mol['O2'][i]  + (F_in['O2'][i]  - n_mol['O2'][i]/tau_head[i]  - r_SOB[i]*Vgas[i]/PM['O2']/R_C)*dt
        n_mol['SX'][i+1]  = n_mol['SX'][i]  + (F_in['SX'][i]  - n_mol['SX'][i]/tau_head[i]  + r_SOB[i]*Vgas[i]/PM['SX'])*dt
        n_mol['H2O'][i+1] = n_mol['H2O'][i] + (F_in['H2O'][i] - n_mol['H2O'][i]/tau_head[i] + r_SOB[i]*Vgas[i]/PM['SX'])*dt
        
        F_out['H2S'][i+1] = n_mol['H2S'][i+1]/tau_head[i+1]
        F_out['O2'][i+1]  = n_mol['O2'][i+1]/tau_head[i+1]    
        F_out['H2O'][i+1] = n_mol['H2O'][i+1]/tau_head[i+1]
        F_out['CO2'][i+1] = F_in['CO2'][i+1]
        F_out['CH4'][i+1] = F_in['CH4'][i+1]
        F_out['N2'][i+1]  = F_in['N2'][i+1]
    
        F_out['total'][i+1] = F_out['H2S'][i+1] + F_out['O2'][i+1] + F_out['H2O'][i+1] + F_out['CO2'][i+1] + F_out['CH4'][i+1] + F_out['N2'][i+1]
        
        x_out['H2S'][i+1]  = F_out['H2S'][i+1]*tau_head[i+1] / (F_out['H2O'][i+1]*tau_head[i+1]  + F_out['H2S'][i+1]*tau_head[i+1] + F_out['O2'][i+1] *tau_head[i+1]+  F_out['CH4'][i+1]*tau_head[i+1] +  F_out['CO2'][i+1]*tau_head[i+1])
        if x_out['H2S'][i+1]*1e+6 > H2S_treshold: 
            if i == 0:
                print('At initial condition H2S concentration is above thresold. \n Unfortunately, even if we can work on it, it is still not possible to change the past. \n Sorry not sorry, we move on.')
                break
            print('*** CRITICAL H2S CONCENTRATION *** \n t = ', t[i+1], 'h \t H2S = ', x_out['H2S'][i+1]*1e+6, 'ppm \t O2,in = ', F_in['O2'][i-1], 'mol/h')
            F_in['O2'][i-1] = min(F_in_max,(F_in['O2'][i-1]*1.2))# Increase O2 influent to have a larger concentration of O2 in the headspace at t=i -> r_SOB[i] can  be increased
            n_mol['O2'][i]  = n_mol['O2'][i-1]  + (F_in['O2'][i-1]  - n_mol['O2'][i-1] /tau_head[i-1]  - r_SOB[i-1]*Vgas[i-1]/PM['O2']/R_C)*dt
            F_out['O2'][i]  = n_mol['O2'][i]/tau_head[i] 
            if F_in['O2'][i-1] >= F_in_max:
                print('Maximum O2 influent reached. Need to anticipate.')
                j = i-2
                while x_out['H2S'][i+1]*1e+6 >= H2S_treshold:                
                    for zz in range(j,i+1): # We need to anticipate the influent of O2 to have a larger concentration of O2 in the headspace at t=i -> iteration up to range(t=i+1) so that stops at t=i)
                        F_in['O2'][zz] = F_in_max
                        
                        x_out['H2S'][zz]  = F_out['H2S'][zz]*tau_head[zz] / (F_out['H2O'][zz]*tau_head[zz]  + F_out['H2S'][zz]*tau_head[zz] + F_out['O2'][zz] *tau_head[zz]+  F_out['CH4'][zz]*tau_head[zz] +  F_out['CO2'][zz]*tau_head[zz] + F_out['N2'][zz]*tau_head[zz])
                        x_out['O2'] [zz]  = F_out['O2'][zz]*tau_head[zz] /  (F_out['H2O'][zz]*tau_head[zz]  + F_out['H2S'][zz]*tau_head[zz] + F_out['O2'][zz] *tau_head[zz]+  F_out['CH4'][zz]*tau_head[zz] +  F_out['CO2'][zz]*tau_head[zz] + F_out['N2'][zz]*tau_head[zz])
                        
                        c_out['H2S'][zz] = x_out['H2S'][zz] * P/R/T # mol/L
                        c_out['O2'][zz]  = x_out['O2'][zz]  * P/R/T # mol/L
                        
                        w_out['H2S'][zz] = c_out['H2S'][zz] * PM['H2S'] * 1000 # g/m3
                        w_out['O2'][zz]  = c_out['O2'][zz]  * PM['O2'] * 1000 # g/m3
                        
                        r_SOB[zz] = (k * w_out['H2S'][zz]**alfa * w_out['O2'][zz]**beta) # gS/m3/h = mgS/L/h

                        n_mol['H2S'][zz+1] = n_mol['H2S'][zz] + (F_in['H2S'][zz] - n_mol['H2S'][zz]/tau_head[zz] - r_SOB[zz]*Vgas[zz]/PM['SX'])*dt
                        n_mol['O2'][zz+1]  = n_mol['O2'][zz]  + (F_in['O2'][zz]  - n_mol['O2'][zz]/tau_head[zz]  - r_SOB[zz]*Vgas[zz]/PM['O2']/R_C)*dt
                        n_mol['SX'][zz+1]  = n_mol['SX'][zz]  + (F_in['SX'][zz]  - n_mol['SX'][zz]/tau_head[zz]  + r_SOB[zz]*Vgas[zz]/PM['SX'])*dt
                        n_mol['H2O'][zz+1] = n_mol['H2O'][zz] + (F_in['H2O'][zz] - n_mol['H2O'][zz]/tau_head[zz] + r_SOB[zz]*Vgas[zz]/PM['SX'])*dt
                        
                        F_out['H2S'][zz+1] = n_mol['H2S'][zz+1]/tau_head[zz+1]
                        F_out['O2'][zz+1]  = n_mol['O2'][zz+1]/tau_head[zz+1]    
                        F_out['H2O'][zz+1] = n_mol['H2O'][zz+1]/tau_head[zz+1]
                        F_out['CO2'][zz+1] = F_in['CO2'][zz+1]
                        F_out['CH4'][zz+1] = F_in['CH4'][zz+1]
                        F_out['N2'][zz+1]  = F_in['N2'][zz+1]
                        F_out['total'][i+1] = F_out['H2S'][i+1] + F_out['O2'][i+1] + F_out['H2O'][i+1] + F_out['CO2'][i+1] + F_out['CH4'][i+1] + F_out['N2'][i+1]
        
                        x_out['H2S'][i+1]  = F_out['H2S'][i+1]*tau_head[i+1] / (F_out['H2O'][i+1]*tau_head[i+1]  + F_out['H2S'][i+1]*tau_head[i+1] + F_out['O2'][i+1] *tau_head[i+1]+  F_out['CH4'][i+1]*tau_head[i+1] +  F_out['CO2'][i+1]*tau_head[i+1])
                        print('Actual time is t = ', t[i], 'h')
                        print('t = ', t[zz])
                        if zz == i and x_out['H2S'][i+1]*1e+6 < H2S_treshold: # If cycle activates only at the last iteration of the loop
                            print('H2S actual is ', x_out['H2S'][i+1]*1e+6, 'ppm -> OK, move on')
                            break
                        elif zz == i and x_out['H2S'][i+1]*1e+6 >= H2S_treshold:
                            print('H2S actual is ', x_out['H2S'][i+1]*1e+6, 'ppm -> still too high, increase the O2 flow earlier')
                            j = j -1
                            if j == 0:
                                print('O2 flow is already at the maximum since the beginning. \n Again, we cannot change the past.')
                                break                        
                break
        else:
            break
                        
print('Simulation ended')
print('time for procedure: ', time.time() - t0, 'sec')

digester_out = pd.DataFrame()
digester_out['t'] = t
digester_out['t_out'] = t+ tau_head
digester_out['F_CH4'] = F_out['CH4']         # [mol/h]   - Methane flow                 
digester_out['F_CO2'] = F_out['CO2']         # [mol/h]   - CO2 flow 
digester_out['F_H2S'] = F_out['H2S']         # [mol/h]   - H2S flow 
digester_out['F_H2O'] = F_out['H2O']         # [mol/h]   - Water flow 
digester_out['F_O2']  = F_out['O2']          # [mol/h]   - Oxygen flow 
digester_out['F_N2']  = F_out['N2']          # [mol/h]   - Nitrogen flow
digester_out['F_total'] = digester_out['F_CH4'] + digester_out['F_CO2'] + digester_out['F_H2S'] + digester_out['F_H2O'] + digester_out['F_O2'] + digester_out['F_N2'] # [mol/h]   - Total flow
digester_out['x_CH4'] = F_out['CH4']/digester_out['F_total']    #[mol/mol] - Outlet molar fraction
digester_out['x_CO2'] = F_out['CO2']/digester_out['F_total']    #[mol/mol] - Outlet molar fraction
digester_out['x_H2S'] = F_out['H2S']/digester_out['F_total']    #[mol/mol] - Outlet molar fraction
digester_out['x_H2O'] = F_out['H2O']/digester_out['F_total']    #[mol/mol] - Outlet molar fraction
digester_out['x_O2']  = F_out['O2']/digester_out['F_total']     #[mol/mol] - Outlet molar fraction
digester_out['x_N2']  = F_out['N2']/digester_out['F_total']     #[mol/mol] - Outlet molar fraction
digester_out['n_SX']  = F_out['SX']                    #[mol]     - Sulfur content
digester_out['F_CH4_in'] = F_in['CH4']         # [mol/h]   - Methane flow
digester_out['F_CO2_in'] = F_in['CO2']         # [mol/h]   - CO2 flow
digester_out['F_H2S_in'] = F_in['H2S']         # [mol/h]   - H2S flow
digester_out['F_H2O_in'] = F_in['H2O']         # [mol/h]   - Water flow
digester_out['F_O2_in']  = F_in['O2']          # [mol/h]   - Oxygen flow
digester_out['F_N2_in']  = F_in['N2']          # [mol/h]   - Nitrogen flow
digester_out['F_total_in'] = digester_out['F_CH4_in'] + digester_out['F_CO2_in'] + digester_out['F_H2S_in'] + digester_out['F_H2O_in'] + digester_out['F_O2_in'] + digester_out['F_N2_in'] # [mol/h]   - Total flow
digester_out['x_CH4_in'] = F_in['CH4']/digester_out['F_total_in']    #[mol/mol] - Outlet molar fraction
digester_out['x_CO2_in'] = F_in['CO2']/digester_out['F_total_in']    #[mol/mol] - Outlet molar fraction
digester_out['x_H2S_in'] = F_in['H2S']/digester_out['F_total_in']    #[mol/mol] - Outlet molar fraction
digester_out['x_H2O_in'] = F_in['H2O']/digester_out['F_total_in']    #[mol/mol] - Outlet molar fraction
digester_out['x_O2_in']  = F_in['O2']/digester_out['F_total_in']     #[mol/mol] - Outlet molar fraction
digester_out['x_N2_in']  = F_in['N2']/digester_out['F_total_in']     #[mol/mol] - Outlet molar fraction
digester_out['r_SOB'] = r_SOB              # [mg/L/h] - Reaction rate
digester_out['r_SOB'] = digester_out['r_SOB'].fillna(method='ffill') # [mg/L/h] - Reaction rate
digester_out['efficiency'] =(1-digester_out['F_H2S']/digester_out['F_H2S_in'])*100 # [-] - Efficiency

print(digester_out)
t = digester_out['t_out']

plt.figure()
plt.subplots_adjust(hspace=0.5)
plt.subplot(2,1,1)
plt.suptitle('Headspace behavior')
plt.title('$H_2S$ Flowrate')
plt.plot(t, F_in['H2S'], '--', label = 'in',  linewidth=2)
plt.plot(t, F_out['H2S'], label = 'out', linewidth=2)

plt.grid(True); plt.xlabel('hours'); plt.ylabel('mol/h'); plt.legend(loc='best')
plt.subplot(2,1,2)
plt.title('$O_2$ Flowrate')
plt.plot(t, F_in_max*np.ones(len(t)), '-.', label = 'Maximum Limit', color = 'red',  linewidth=2.5)
plt.plot(t, F_in['O2'], '--', label = 'in',  linewidth=2)
plt.plot(t, F_out['O2'], label = 'out', linewidth=2)
plt.grid(True); plt.xlabel('Time [h]'); plt.ylabel('mol/h'); plt.legend(loc='best')

plt.figure()
plt.suptitle('Biogas Composition')
plt.subplot(3,1,1)
plt.plot(t, digester_out['x_CH4']*100, label = '$CH_4$', linewidth=2)
plt.plot(t, digester_out['x_CO2']*100, label = '$CO_2$', linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('Molar Fraction [%]'); plt.legend()
plt.subplot(3,1,2)
plt.plot(t, digester_out['x_H2S_in']*1e+6, '--', label = '$x_{H_2S,in}$', linewidth=2)
plt.plot(t, digester_out['x_H2S']*1e+6, label = '$x_{H_2S,out}$', linewidth=2)
plt.plot(t, np.ones(len(t))*H2S_treshold, '--', color = 'red', label = 'Threshold', linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('Molar Fraction [ppm]'); plt.legend()
plt.subplot(3,1,3)
plt.plot(t, digester_out['x_O2']*100, label = '$O_2$', linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('Molar Fraction [%]'); plt.legend()
plt.legend()

plt.figure()
plt.title('Microaeration Effect')
plt.plot(t, digester_out['x_H2S']*1e+6, label = '$x_{H_2S,out}$', linewidth=2)
plt.plot(t, np.ones(len(t))*H2S_treshold, '--', color = 'red', label = 'Threshold', linewidth=2)
plt.grid(which='both'); plt.xlabel('hours'); plt.ylabel('$H_2S$ Mole Fraction [ppm]'); plt.legend()
plt.show()

