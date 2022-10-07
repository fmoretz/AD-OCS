import pandas as pd


simname = "provaADM1"
simname = input("Insert the name of the simulation:")
print("Data from:",simname)
folder  =  r'C:\Users\fmoretta\OneDrive - Politecnico di Milano\SUPER\PhD\Progetti Magistrali\Federico Rocca Thesis\Tecnologia\Data\Worksheets'
folder  = r'C:\Users\fede1\OneDrive - Politecnico di Milano\Federico Rocca Thesis\Tecnologia\Data\Worksheets'
reading_path =  folder + "\\" + simname + "py"+ ".xlsx"

colnames = ["HRT","S1","XT", "S2", "X1", "X2", "Z", "C","CO2","B", "pH", "q_C", "P_C", "q_CH4"]

T1 = pd.read_excel(reading_path, sheet_name = "SS_Values",header = None, names = colnames, skiprows = 1)
T2 = pd.read_excel(reading_path, sheet_name = "Influent", header = 0)
