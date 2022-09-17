# Reactor Conditions
import numpy as np
from PhysConstants import*

D    = 0.05    # [1/d] Dilution rate
T    = 35      # [°C]
Pt   = 1       # [atm]
alfa = 1
V_liq = 3400    # [m3] Volume of the reactor
V_gas = 300     # [m3] Volume of the headspace
V_ad = V_liq + V_gas

C_d = [0.1, 0.1]

sol  = np.exp(-159.854 + 8741.68/(T+273.15) + 21.6694*np.log(T+273.15) - 0.00110261*(T+273.15))   # [-] mole fraction of dissolved CO2 in water at [T]°C
KH = 1/sol                                                        # [atm]        - Henry's constant CO2 at [T]°C - Partial Pressure Relation
KH = 1/(KH/55342)                                                 # [mmol/L/atm] - Henry's constant CO2 at [T]°C - Concentrations Relation


# H2S Data
H_S  = 0.1367719*T + 2.0180702               # [atm]        - Henry's constant CO2 at [T]°C - Partial Pressure Relation: Linear regression on Perry's data
H_S  = 0.00030308*(T**2) + 0.11276471*T +2.44557423
KH_S = 1/(H_S*100/55342)                    # [mmol/L/atm] - Henry's constant CO2 at [T]°C - Concentrations Relation
H_S  = 1/KH_S                               # [atm L/mmol] 
