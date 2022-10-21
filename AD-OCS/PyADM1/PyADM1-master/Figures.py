import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


influent_state = pd.read_csv('digester_influent.csv')
effluent_state = pd.read_csv("dynamic_out.csv")
gas = pd.read_csv("gas_flows.csv")
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

plt.figure(5)
plt.subplot(2,1,1)
plt.plot(timespan, gas['q_ch4'], label='CH4')
# plt.plot(timespan, gas['q_co2'], label='CO2')
#plt.plot(timespan, gas['q_h2'], label='H2')
plt.legend(loc='best')
plt.subplot(2,1,2)
plt.plot(timespan, gas['q_co2'], label='CO2')
plt.legend(loc='best')
plt.show()