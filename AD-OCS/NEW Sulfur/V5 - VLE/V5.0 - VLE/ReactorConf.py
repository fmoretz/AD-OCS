# Reactor Conditions
import numpy as np
from PhysConstants import*

D    = 0.05     # [1/d] Dilution rate
T    = 35       # [°C]
Pt   = 1        # [atm]
alfa = 1
V_liq = 3400    # [m3] Volume of the reactor
V_gas = 300     # [m3] Volume of the headspace
V_ad = V_liq + V_gas

C_d = [0.1, 0.1]

P_dig = 1 # [atm] Pressure of the digester

