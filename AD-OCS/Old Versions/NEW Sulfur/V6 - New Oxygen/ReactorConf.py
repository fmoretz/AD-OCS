# Reactor Conditions
import numpy as np
from PhysConstants import*

D    = 0.05     # [1/d] Dilution rate
T    = 35       # [°C]
Pt   = 1        # [atm]
alfa = 1

C_d = [0.1, 0.1]

P_dig = 1 # [atm] Pressure of the digester

Dr = 20 # [m] Diameter of the digester
H0 = 10 # [m] Initial height of liquid in the digester
SR = np.pi/4 * Dr**2 # m2    - reactor surface
hmax = 14   # [m] Maximum height of liquid in the digester
hmin = 5    # [m] Minimum height of liquid in the digester

V_reactor = 3400    # [m3] Volume of the reactor
V_headspace = 450     # [m3] Volume of the headspace
V_ad = V_reactor + V_headspace