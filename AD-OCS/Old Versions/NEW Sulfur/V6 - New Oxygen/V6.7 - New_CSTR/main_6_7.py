''' Current version of AD_OCS model. Define the input according to the defined unit of measure from the spreadsheet.
    Define the timespan in hours.'''

import math
from pprint import pprint as pp
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

from SS_Algebraic import*
from dataimport import*
from PhysConstants  import*

from functions import gompertz, growth_SRB, f_deviations, deviations_check, level_t, headspace_dynamics_discr, logistic_deviations,  headspace_dynamics_discr_SS
from models import AD_OCS, AMOCO_HN
from mix_real_gas import f_VL_realgas


# System definition
d_start = 0        # [h] - Start time
d_end   = 60       # [h] - End time
hours   = 0.2     # [h] - Discretization time
X_CRIT = 150 # ppm

n_times = int((d_end-d_start)/hours)+1  # Number of time steps

print('***Intervals of {hours} hours. \n {n_times} time steps***'.format(hours=hours, n_times=n_times))

t_span = np.linspace(d_start,d_end,n_times) # time span
t_span_d = t_span/24 # time span in days

y_influent_changes = f_deviations(t_span_d, T3.index.values, y_in_0) # Get the deviated influent values at each timestamp
y_influent = y_influent_changes[:,0:5] # Remove the Q_in values

Q_in = y_influent_changes[:,5] # Remove the y_in values

# --------------------------------------------------------------------------------------------
# Ode Integration of AMOCO_HN: use to get X2 
# --------------------------------------------------------------------------------------------

y0= [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]] # initial conditions esxkcdlished from SS

YOUT_pre = odeint(AMOCO_HN, y0, t_span_d, hmax = (t_span_d[1]-t_span_d[0]), args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values))
# X2_pre = YOUT_pre[:,2]
S2_pre = YOUT_pre[:,5]
print('************** AMOCOHN OK *******************')

# --------------------------------------------------------------------------------------------
# Ode Integration of the AD_OCS_Model: uses the previous X2 to get rho
# --------------------------------------------------------------------------------------------

y0 = [SSTATE[0], SSTATE[1], SSTATE[2], SSTATE[3], SSTATE[4], SSTATE[5], SSTATE[6]] # initial conditions esxkcdlished from SS

YOUT= odeint(AD_OCS, y0, t_span_d, hmax = (t_span_d[1]-t_span_d[0]), args=(alfa, mu_max, Ks, KI2, KH, Pt, kLa, D, k, kd, N_bac, N_S1, y0[2], t_span[0], y_in_0, T3.index.values, t_span_d))

# Get results
XT = YOUT[:,0]              # [gCOD/L] - Particulate 
X1 = YOUT[:,1]              # [g/L]    - Acidogenics  Bacteria  
X2 = YOUT[:,2]              # [g/L]    - Methanogenic Bacteria
Z  = YOUT[:,3]              # [mmol/L] - Total Alkalinity
S1 = YOUT[:,4]              # [g/L]    - Organic Soluble Substrate
S2 = YOUT[:,5]              # [mmol/L] - VFA dissolved
C  = YOUT[:,6]              # [mmol/L] - Inorganic Carbon Dissolved

print('************** AD_OCS OK *******************')

# Solver Output: from all the variables from the ones of the ODE
mu1 = np.empty(len(XT))
mu2 = np.empty(len(XT))
CO2 = np.empty(len(XT))
B   = np.empty(len(XT))
phi = np.empty(len(XT))
p_C = np.empty(len(XT))
q_C = np.empty(len(XT))
q_M_u = np.empty(len(XT))
pH  = np.empty(len(XT))

for x in range(len(t_span)):
    mu1[x] = mu_max[0]*(S1[x]/(S1[x]+Ks[0]))                     # [1/d]      - Specific Growth Rate for X1 (Monod)
    mu2[x] = mu_max[1]*(S2[x]/(S2[x]+Ks[1]+S2[x]**2/KI2))        # [1/d]      - Specific Growth Rate for X2 (Haldane)
    CO2[x] = C[x] + S2[x] - Z[x]                                 # [mmol/L]   - Dissolved CO2
    B[x]   = Z[x] - S2[x]                                        # [mmol/L]   - Alkalinity
    phi[x] = CO2[x] + KH*Pt + k[5]/kLa*mu2[x]*X2[x]
    p_C[x]  = (phi[x] - (phi[x]**2- 4*KH*Pt*CO2[x])**0.5)/(2*KH) # [atm]      - CO2 Partial Pressure
    q_C[x] = kLa*(CO2[x] - KH*p_C[x])                            # [mmol/L/d] - CO2 Outlet Molar Flow
    q_M_u[x] = k[5]*mu2[x]*X2[x]                                 # [mmol/L/d] - CH4 Outlet Molar Flow
    pH[x]  = np.real(-np.log10(Kb*CO2[x]/B[x]))                  # [-]        - System pH


# Sulfur and Oxygen Influence Evaluation

Xs     = np.zeros(len(t_span))                                   # [g/L]      - Sulfate Reducing Bacteria
Ss     = np.zeros(len(t_span))                                   # [mmol/L]      - Sulfur dissolved
y_S    = np.zeros(len(t_span))                                   # [-]        - Sulfur Mole fraction in gas phase
Ss_max = np.zeros(len(t_span))                                   # [mmol/L]      - Maximum Sulfur dissolved concentration
Xs_max = np.zeros(len(t_span))                                   # [g/L]      - Maximum Sulfate Reducing Bacteria concentration (Gompertz Asymptote)
mu_srb = np.zeros(len(t_span))                                   # [g/L/d]    - Gompertz parameter for SRB growth
growth_rate = np.zeros(len(t_span))                              # [g/L/d]    - Gompertz derivative: states the growth rate of SRB

dXsdt = np.zeros([len(t_span),len(t_span)])                     # [g/L/d]    - SRB growth rate matrix, preallocation
lam = 0                                                         # [-]        - Lag Time as defined by Gompertz

index = 0
for j in range(len(t_span)):
    # Iterate over each time step
        
    Ss_max[j] = frac_sulfur*1000/64*S2[j]                   #*y_influent[j,4]/y_influent[0,4]
    Xs_max[j] = Y_srb/(1-Y_srb)*Ss_max[j]     
    
    # Gompertz function for microbial population and dissolved sulfur
    mu_srb[j] = (- X2[0] + X2[j])/max(1e-14,t_span[j]-t_span[0])      # [g/L/d]    - Gompertz parameter for SRB growth
    Xs[j]  = gompertz(t_span[j], Xs_max[j], mu_srb[j], lam)           # [g/L]      - Sulfate Reducing Bacteria - Gompertz
    Ss[j]  = Xs[j]*(1-Y_srb)/(Y_srb)                                  # [mmol/L]      - Sulfur dissolved concentration
    
    growth_rate[j] = mu2[j]*Xs[j]/X2[j]                       # [g/L/d]    - Get the growth rate of SRB at each time step as the maximum of the possible rates

    y_S[j]   = (KH_S*Ss[j])/P_dig                             # [-]        - Sulfur Mole fraction in gas phase (G/L equilibrium)


I_Ss = 1-KI_SS*Ss                                             # [-]   - Sulfate Ion concentration  
q_M  = q_M_u*I_Ss
q_tot = q_C + q_M                                             # [mmol/L/d] - Outlet global molar flow  
q_S = y_S*q_tot                                               # [mmol/L/d] - Sulfur Outlet Specific Molar Flow

### Assess the liquid volume dynamics ###
# Smooth the step change of the liquid influent
Q_in_real = np.ones(len(Q_in))*Q_in[0]
loc_index = 0

for i in range(len(Q_in)):
    t=t_span[i]    
    if t < T3.index.values[0]*24:
        Q_in_real[i] = Q_in[i]
    elif t >= T3.index.values[-1]*24:
        Q_in_0_loc = Q_in_real[t_span<T3.index.values[-1]*24][-1]
        Q_in_real[i] = logistic_deviations(t, T3.Qin.iloc[-1]*Q_in[0], 1, T3.index.values[-1]*24, Q_in_0_loc)
    else:
        Q_in_0_loc = Q_in_real[t_span<T3.index.values[loc_index]*24][-1]
        Q_in_real[i] = logistic_deviations(t, T3.Qin.iloc[loc_index]*Q_in[0], 1, T3.index.values[loc_index]*24, Q_in_0_loc)
        if t >= T3.index.values[loc_index+1]*24:
            loc_index += 1
Q_df = pd.DataFrame(index=t_span)
Q_df = pd.DataFrame()
Q_df["Qin"] = Q_in
Q_df["Qin_real"] = Q_in_real
Q_df["Qin_smooth"] = Q_df["Qin_real"].rolling(window=10).mean()
Q_df["Qin_smooth"].fillna(method='bfill', inplace=True)
Q_in_real = Q_df["Qin_smooth"].values

h = np.zeros(len(t_span))
h0 = h_SS
t_change = t_span[0]
index = 0
V_liq = np.zeros(len(t_span))                                 # [m3] - Liquid Volume
contr = 0
for i in range(len(t_span)):
    t = t_span[i]  
    if t < T3.index.values[0]*24:
        h[i] = h0
    elif t >= T3.index.values[-1]*24:
        if contr == 0:
            h0_loc = h[t_span<T3.index.values[index]*24][-1]
            h0_loc = level_t(t, D, Q_in_real[i], SR, h0_loc, t_change)  
            contr = 1
        t_change = T3.index.values[-1]*24
        h[i] = level_t(t, D, Q_in_real[i], SR, h0_loc, t_change)     
    else:    
        h0_loc = h[t_span<T3.index.values[index]*24][-1]
        t_change = T3.index.values[index]*24
        h[i] = level_t(t, D, Q_in_real[i], SR, h0_loc, t_change)     
        if t >= T3.index.values[index+1]*24:
            index += 1               
    
    if h[i] > hmax:
        print('!!! Level is too high !!!')
        input('Press Enter to continue')
    elif h[i] < hmin:
        print('!!! Level is too low !!!')
        input('Press Enter to continue')
h_df = pd.DataFrame(index=t_span_d)
h_df = pd.DataFrame()
h_df["h"] = h
h_df["h_smooth"] = h_df["h"].rolling(window=5).mean()
h_df["h_smooth"].fillna(method='bfill', inplace=True)
h = h_df["h_smooth"].values
V_liq = np.pi*Dr**2*h/4                             # [m3] - Liquid Volume

#### Evaluation of the G/L equilibrium effects on the system ####

n_species = 4                                                 # [-]       - Number of species in the system

F_i = np.zeros([len(t_span),n_species])                       # [-] - Outlet Molar Fraction - In of the flash
F_M = np.zeros(len(t_span))                                   # [-] - Methane Flow
F_C = np.zeros(len(t_span))                                   # [-] - Carbon Flow
F_S = np.zeros(len(t_span))                                   # [-] - Sulfur Flow
F_W = np.zeros(len(t_span))                                   # [-] - Water Flow

for i in range(len(t_span)):  
    F_W[i] = water_percentage*Q_in_real[i]*rho_W/18*1000/24             # [mol/h]   - Water Flow 
    F_M[i] = q_M[i]*V_liq[i]/24                                    # [mol/h]   - Methane Flow 
    F_C[i] = q_C[i]*V_liq[i]/24                                    # [mol/h]   - CO2 Flow
    F_S[i] = q_S[i]*V_liq[i]/24                                    # [mol/h]   - Sulfur Flow

    F_i[i] = np.array([F_M[i], F_C[i], F_S[i], F_W[i]])         # [mol/h] - Outlet Molar Flow - In of the flash

z_i = np.zeros([len(t_span),n_species])                         # [-] - Outlet Molar Fraction - In of the flash
x_i = np.zeros([len(t_span),n_species])                         # [-] - Outlet Molar Fraction - Liquid phase of the flash
y_i = np.zeros([len(t_span),n_species])                         # [-] - Outlet Molar Fraction - Vapor phase of the flash

K = np.zeros([len(t_span),n_species])                         # [-] - Equilibrium Constant
num = np.zeros([len(t_span),n_species])                       # [-] - Numerator of the equilibrium equation
den = np.zeros([len(t_span),n_species])                       # [-] - Denominator of the equilibrium equation

N_L_tot = np.zeros(len(t_span))                               # [mol/h] - Total Liquid Molar Flow
N_V_tot = np.zeros(len(t_span))                               # [mol/h] - Total Vapor Molar Flow

N_L = np.zeros([len(t_span),n_species])                       # [mol/h] - Liquid Phase Molar Flow
N_V = np.zeros([len(t_span),n_species])                       # [mol/h] - Vapor Phase Molar Flow
Q_L = np.zeros([len(t_span),n_species])                       # [m3/h] - Liquid Phase volumetric flow
Q_V = np.zeros([len(t_span),n_species])                       # [m3/h] - Vapor Phase volumetric flow

alpha_flash = np.zeros(len(t_span))                           # [-]     - 
MW_vett = [MW_M, MW_C, MW_S, MW_W]                            # [g/mol] - Molecular Weight Vector
rho_vett = [rho_M, rho_C, rho_S, rho_W]                       # [g/L]   - Density Vector
law = ('H', 'H', 'H', 'R')                                    # [-]     - Define the law to be used for the equilibrium calculation (R for Raoult's Law, H for Henry's Law)
for t in range(len(t_span)):  
    F_in = sum(F_i[t])                                        # [mol/h] - Total Inlet Molar Flow
    for i in range(n_species):        
        z_i[t,i] = F_i[t,i]/F_in                              # [-]     - Outlet Molar Fraction - In of the flash        
      
    VLsolution = f_VL_realgas(z_i[t], T, P_dig)               # [-]     - Solve the G/L equilibrium
    alpha_flash[t] = VLsolution[0]                            # [-] - Vapor fraction      
    K[t]= VLsolution[1]         
    
    N_V_tot[t] = F_in*alpha_flash[t]                                # [mol/h] - Vapor Molar Flow
    N_L_tot[t] = F_in - N_V_tot[t]                                  # [mol/h] - Liquid Molar Flow
    for i in range(n_species):
        x_i[t,i] = z_i[t,i]/(1+alpha_flash[t]*(K[t,i]-1))          # [-]     - Outlet Molar Fraction - Liquid phase of the flash
        y_i[t,i] = K[t,i]*x_i[t,i]                                 # [-]     - Outlet Molar Fraction - Vapor phase of the flash
        N_L[t,i] = N_L_tot[t]*x_i[t,i]                             # [mol/h] - Liquid Molar Flow of species i
        N_V[t,i] = N_V_tot[t]*y_i[t,i]                             # [mol/h] - Vapor Molar Flow of species i
        Q_L[t,i] = N_L[t,i]*MW_vett[i]/rho_vett[i]/1000            # [m3/h]  - Liquid Volumetric Flow of species i
        Q_V[t,i] = N_V[t,i]*(Rgas_m3_atm_K*(T+273.15)/P_norm)           # [m3/h]  - Vapor Volumetric Flow of species i


#### Evaluation of the Oxygen effect equilibrium effects on the system: CSTR behaviour ####
# IN: N_V[M, C, S, W] + N_in_O2
# OUT: N_out[M, C, S, W, O]


RC = 1.8                                                         # [gS/gO2]     - Stoichiometry of the reaction
alpha = 1                                                       # [-] - Power for sulfur concentration
beta = 0.3                                                     # [-] - Power for oxygen concentration
k_SOB = 0.8                                                    # [1/m3/h*(gS*gO)^proper exp] - Reaction rate constant
r_vett = np.zeros(len(t_span))                                  # [mg/L/min] - Reaction rate
S_vett = np.zeros((len(t_span),4))                              # [mg/L] - Sulfur concentration

V_gas = ((V_reactor - V_liq) + V_headspace)

tau_headspace = V_gas/Q_V.sum(1)                            # [h] - Time to fill the headspace

Q_in = 1                                                   # [m3/h] - Inlet volumetric flow - Gas (air or oxygen) defined in the input section
if air == 1:
    N_in_O2 = np.ones(len(t_span))*Q_in*0.21*P_norm/(Rgas_m3_atm_K*(T_norm+273.15))    # [mol/h] - Inlet molar flow - Oxygen
    N_in_N2 = np.ones(len(t_span))*Q_in*0.79*P_norm/(Rgas_m3_atm_K*(T_norm+273.15))    # [mol/h] - Inlet molar flow - Nitrogen
elif air == 2:
    N_in_O2 = np.ones(len(t_span))*Q_in*P_norm/(Rgas_m3_atm_K*(T_norm+273.15))         # [mol/h] - Inlet molar flow - Oxygen
    N_in_N2 = np.zeros(len(t_span))
else:
    print('Error: air must be 1 or 2')

N_in_gas = {'O2': N_in_O2,
            'N2': N_in_N2}# [mol/h] - Oxygen Flow with respect to sulfur flow ==========================================OXYGEN FLOW

F_OUT_0 = {'H2S': N_V[0,2]/2,
            'CH4': N_V[0,0],
            'CO2': N_V[0,1],
            'O2':  N_in_O2[0]/2,
            'H2O': N_V[0,3],
            'SX':  0,
            'N2':  N_in_N2[0]
            } # [mol/h] - Initial flow, SX is [mol]

F_out_SS, n_SX_SS = headspace_dynamics_discr_SS(np.linspace(0,100,1000), N_V, N_in_gas, F_OUT_0, V_gas, tau_headspace, RC, alpha, beta, k_SOB, P_dig, T, T_norm, P_norm)
F_OUT_SS = {k: v[-1] for k, v in F_out_SS.items()}
F_OUT_SS['SX'] = n_SX_SS[-1]

F_out, F_in, r_SOB = headspace_dynamics_discr(t_span, N_V, N_in_gas, F_OUT_SS, V_gas, tau_headspace, RC, alpha, beta, k_SOB, P_dig, T, T_norm, P_norm)

digester_out = pd.DataFrame()
digester_out['t'] = t_span
digester_out['t_out'] = t_span + tau_headspace
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


digester_out.head()
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
plt.grid(True); plt.xlabel('hours'); plt.ylabel('Molar Fraction [ppm]'); plt.legend()
plt.subplot(3,1,3)
plt.plot(t, digester_out['x_O2']*100, label = '$O_2$', linewidth=2)
plt.grid(True); plt.xlabel('hours'); plt.ylabel('Molar Fraction [%]'); plt.legend()
plt.legend()

plt.show()