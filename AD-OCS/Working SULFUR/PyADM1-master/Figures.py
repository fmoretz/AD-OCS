import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

influent_state = pd.read_csv("digester_influent.csv")
effluent_state = pd.read_csv("dynamic_out.csv")
timespan = influent_state['time']
particu = effluent_state['X_xc']
S_SS_su = effluent_state['S_su']

plt.figure(1)
plt.plot(timespan, particu, label='Particulate')
plt.figure(2)
plt.plot(timespan, effluent_state['S_su'], label='S_su')
plt.plot(timespan, effluent_state['X_su'], label='X_su')
plt.legend(loc='best')
plt.figure(3)
plt.plot(timespan, effluent_state['S_aa'], label='S_aa')
plt.plot(timespan, effluent_state['X_aa'], label='X_aa')
plt.legend(loc='best')
plt.figure(4)
plt.plot(timespan, effluent_state['S_fa'], label='S_fa')
plt.plot(timespan, effluent_state['X_fa'], label='X_fa')
plt.legend(loc='best')
plt.show()