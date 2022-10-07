import pandas as pd

# simname = input("Insert the name of the simulation:")
simname = "AMOCO_HN"
print("Data from:",simname)

folder  =  r'C:\Users\fede1\OneDrive - Politecnico di Milano\Desktop\Tesi\Constant Minimization Rocca'
reading_path =  folder + "/" + simname + "py"+ ".xlsx"

colnames = ["HRT","S1","XT", "S2", "X1", "X2", "Z", "C","CO2","B", "pH", "q_C", "P_C", "q_CH4"]

T1 = pd.read_excel(reading_path, sheet_name = "SS_Values",header = None, names = colnames, skiprows = 1)
T2 = pd.read_excel(reading_path, sheet_name = "Influent", header = 0)

# Get raw data
HRT   = T1["HRT"]
S1    = T1["S1"]  # [g/L] - COD
XT    = T1["XT"]
S2    = T1["S2"]  # [mmol/L] - VFA
X_1   = T1["X1"]
X_2   = T1["X2"]
Z     = T1["Z"]
C     = T1["C"]
CO2   = T1["CO2"]
B     = T1["B"]
pH    = T1["pH"]
q_C   = T1["q_C"]  # [mmol/L/d]
P_C   = T1["P_C"]  # [atm]
q_M   = T1["q_CH4"]

D = 1/HRT

S1_in = float(T2["S1in"])    # [gCOD/L]
S2_in = float(T2["S2in"])    # [mmol/L]
C_in  = float(T2["Cin"])    # [mmol/L]
XT_in = float(T2["XTin"])    # [gCOD/L]
