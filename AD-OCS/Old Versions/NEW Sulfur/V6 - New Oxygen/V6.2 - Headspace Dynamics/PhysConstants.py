# Physical Constants
import numpy as np
from dataimport import *

### Physical Constants ###

Rgas_L_atm_K = 0.08205746 # L*atm/(mol*K)
Rgas_J_mol_K = 8.3144598 # J/(mol*K)
Rgas_m3_atm_K = 0.08205746e-3 # m^3*atm/(mol*K)


MW_APr = 74.08                       # [g/mol] Propionic Acid
MW_AAc = 60.052                      # [g/mol] Acetic Acid
MW_ABu = 88.11                       # [g/mol] Butyric Acid
MW_APe = 102.13                      # [g/mol] Pentanoic Acid
MW_C  = 44                           # [g/mol] CO2
MW_M  = 16                           # [g/mol] CH4
MW_S  = 34                           # [g/mol] H2S
MW_O  = 32                           # [g/mol] O2
MW_W = 18.015                        # [g/mol] H2O

MW_gas  = [MW_M, MW_C, MW_S, MW_O]

rho_W = 997                       # [kg/m^3] Density of water 
rho_M = 0.627                     # [kg/m^3] Density of Methane (engineering toolbox)
rho_C = 1.724                     # [kg/m^3] Density of CO2 (engineering toolbox)
rho_S = 1.505                     # [kg/m^3] Density of H2S (engineering toolbox)

### Biological Parameters ###

N_bac   = 0.08/14*1000*1.55    # Nitrogen content of the biomass [mmolN/gVS]
N_aa    = 0.007*1000*1.55      # Nitrogen content of aminoacids and proteins [mmolN/gVS]
f_aa    = 0.00003              # fraction of aminoacids in the organics substrate S1
f_pr    = 0.62             	   # fraction of proteins in the organic substrate S1
f_xc    = 0.06                 # fraction of composite in the organic substrate S1
f_xc_pr = 0.2                  # fraction of proteins in the composite Xc
N_S1    = (f_aa+f_pr+f_xc*f_xc_pr)*N_aa

Y_srb  = 0.0342*64/1000               # [g/mmol] Yield of SRB population

### Thermodynamic properties ###

sol  = np.exp(-159.854 + 8741.68/(T+273.15) + 21.6694*np.log(T+273.15) - 0.00110261*(T+273.15))     # [-] mole fraction of dissolved CO2 in water at [T]°C (273 K-353 K)
H_C_atm = 1/sol                                                                                     # [atm]        - Henry's constant CO2 at [T]°C - Partial Pressure Relation
KH = 1/(H_C_atm/55342)                                                                            # [mmol/L/atm] - Henry's constant CO2 at [T]°C - Concentrations Relation

# H2S Data

H_S_atm  = (0.1367719*T + 2.0180702)*100         # [atm]        - Henry's constant H2S at [T]°C - Partial Pressure Relation: Linear regression on Perry's data
KH_S = H_S_atm/55342                             # [atm/(mol/m3)] - Henry's constant H2S at [T]°C - Concentrations Relation = [L*atm/mmol]

# Methane Data

sol_M = np.exp(-338.217 + 13281.1/(T+273.15) + 51.9144*np.log(T+273.15) - 0.0425831*(T+273.15))   # [-]   - Mole fraction of dissolved methane in water at [T]°C (273 K-523 K)
H_M_atm = 1/sol_M                                                                                 # [atm] - Henry's constant methane at [T]°C - Partial Pressure Relation

H_atm = [H_M_atm, H_C_atm, H_S_atm, 1]                                                            # [atm] - Henry's constant at [T]°C - Partial Pressure Relation

# Relation for vapor pressures: Perry's Handbook, 8th edition, pages 2-55,2-60 
# P = C1 + C2/T + C3*ln(T) + C4*(T^C5); T [K], P [Pa]
# species [CH4, CO2, H2S]
C1 = [39.205, 140.54, 85.584, 73.649]
C2 = [-1324.4, -4735, -3839.9, -7258.2] 
C3 = [-3.4366, -21.268, -11.199, -7.3037]
C4 = [3.1019E-05, 4.0909E-02, 1.8848E-02, 4.1653E-06]
C5 = [2, 1, 1, 2]

P_sat = np.zeros([len(H_atm)])
for i in range(len(H_atm)):
    P_sat[i] = np.exp(C1[i] + C2[i]/(T+273.15) + C3[i]*np.log(T+273.15) + C4[i]*(T+273.15)**C5[i])  # [Pa] - Vapor pressure of the species at [T]°C
    P_sat[i] = P_sat[i]/101325                                                                      # [atm] - Vapor pressure of the species at [T]°C

H_atm[3] = P_sat[3]                                                                                 # [atm] - Says that for water the Henry's constant is the vapor pressure
