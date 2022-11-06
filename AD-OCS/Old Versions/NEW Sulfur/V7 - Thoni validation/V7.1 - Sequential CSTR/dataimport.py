import pandas as pd
import numpy as np
from pathlib import Path

print('Choose a datasets: \n 1 -> AMOCO_HN \n 2 -> provaADM1 \n 3 -> bsm2 \n 4 -> matlab \n 5 -> thoni')
name_in_0dex = 6 # input("->")

datasets = ["amoco_HN","provaADM1", "bsm2","matlab", "thoni", "thoni_2"]
simname  = datasets[int(name_in_0dex) -1]
print("Data are from:",simname)

filename =simname + '.xlsx'
p = Path.cwd()
reading_path = p.parent / 'Working_data' / filename # Path to the file -  check that it correctly goes back tot the appropriate parent folder from current WD
reading_path = r'C:\Users\fede1\OneDrive - Politecnico di Milano\Documenti\GitHub\AD-OCS\Working_data\thoni_2.xlsx'

colnames = ["HRT", "XT", "S1", "S2", "X1", "X2", "C", "Z", "CO2", "B", "pH", "P_C", "q_C", "q_CH4"]

T1 = pd.read_excel(reading_path, sheet_name = "SS_Values",header = None, names = colnames, skiprows = 1)
T2 = pd.read_excel(reading_path, sheet_name = "Influent", header = 0)
T3 = pd.read_excel(reading_path, sheet_name = "Deviations", header = 0, index_col = 'time')
T4 = pd.read_excel(reading_path, sheet_name = "Reactor", header = 0)

P_norm = 1 # [atm] Normal pressure
T_norm = 20 # [°C] Normal temperature

# Get raw data
HRT   = T1["HRT"]
S1    = T1["S1"]    # [g/L] - COD
XT    = T1["XT"]
S2    = T1["S2"]    # [mmol/L] - VFA
X_1   = T1["X1"]
X_2   = T1["X2"]
Z     = T1["Z"]
C     = T1["C"]
CO2   = T1["CO2"]
B     = T1["B"]
pH    = T1["pH"]
q_C   = T1["q_C"]   # [mmol/L/d]
P_C   = T1["P_C"]   # [atm]
q_M   = T1["q_CH4"]

Dil   = 1/HRT # [1/d] - Vector of dilution rates

# Reactor Configuration 
D    = float(T4["D"])        # [1/d] Dilution rate
T    = float(T4["T"])        # [°C]
Pt   = float(T4["P"])         # [atm]
alfa = float(T4["alpha"]) 

P_dig = Pt # [atm] Pressure of the digester

Dr = float(T4["Dr"])  # [m] Diameter of the digester
H0 = float(T4["H0"])  # [m] Initial height of liquid in the digester
SR = np.pi/4 * Dr**2 # m2    - reactor surface
hmax = float(T4["hmax"])   # [m] Maximum height of liquid in the digester
hmin = float(T4["hmin"])     # [m] Minimum height of liquid in the digester

V_reactor = float(T4["Vr"])     # [m3] Volume of the reactor
V_headspace = float(T4["Vh"])      # [m3] Volume of the headspace
V_ad = V_reactor + V_headspace



# Influent

S1_in_0 = float(T2["S1in"])                              # [g/L]    - COD in solution (soluble COD): S1
S2_in_0 = float(T2["S2in"])                              # [mmol/L] - Volatile Fatty Acids: S2
C_in_0  = float(T2["Cin"] )                              # [mmol/L] - Inorganic Carbon
N_in_0  = float(T2["Nin"] )                              # [mmol/L] - Inorganic Nitrogen
XT_in_0 = float(T2["XTin"])                              # [g/L]    - Particulate COD in solution : XT

Ka   = 1.5e-5       # [?] - Acidity constant of acetic acid
Kb   = 6.5e-7       # [?] - Acidity constant of bicarbonate
Kc   = 4.5e-11      # [?] - Acidity constant of carbonic acid
pH_in_0 = float(T2["pHin"])                              # [-]      - pH
B_in_0  = Kb*C_in_0/(Kb+10**(-pH_in_0))
Z_in_0  = B_in_0 + S2_in_0*Ka/(Ka+10**(-pH_in_0))+N_in_0



Q_in_0 = float(T2["Qin"])                              # [m3/d]  - Influent flowrate
xW_in_0 = float(T2["xW"])                       # [-]    - Fraction standing for percentage of water in the influent
# Sulfur
gammaS_in_0 = float(T2["gammaS"]) 

# Get the array of the influent concentrations and flowrate
y_in_0 = np.array([S1_in_0, S2_in_0, C_in_0, Z_in_0, XT_in_0, Q_in_0, gammaS_in_0, xW_in_0])
# reorder the columns of T3 to get the proper order of the deviations for the functions
T3_m = pd.DataFrame()
T3_m["S1in"] = T3["S1in"]
T3_m["S2in"] = T3["S2in"]
T3_m["Cin"] = T3["Cin"]
T3_m["Nin"] = T3["Nin"]
T3_m["XTin"] = T3["XTin"]
T3_m["Qin"] = T3["Qin"]
T3_m["gammaSin"] = T3["gammaSin"]
T3_m["xWin"] = T3["xWin"]
T3 = T3_m
print(T3)