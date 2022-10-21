import numpy as np
import matplotlib.pyplot as plt

from main_V3 import*

# --------------------------------------------------------------------------------------------
F_M = q_M_new*V_liq # [mol/d] Methane production rate
F_S = q_S_new*V_liq # [mol/d] Sulfur production rate
F_C = q_C_new*V_liq # [mol/d] Carbon dioxide production rate

F = {'CH4': F_M, 'CO2': F_C, 'H2S': F_S, }

# --------------------------------------------------------------------------------------------
Q_M = F_M*T/P_dig*Rgas_m3_atm_K     # [m^3/d] Methane production rate
Q_S = F_S*T/P_dig*Rgas_m3_atm_K     # [m^3/d] Sulfur production rate
Q_C = F_C*T/P_dig*Rgas_m3_atm_K     # [m^3/d] Carbon dioxide production rate

Q = {'CH4': Q_M, 'CO2': Q_C, 'H2S': Q_S, }
Qsum = sum(Q.values())
x_vol = {'CH4': Q_M/(sum(Q.values())), 'CO2': Q_C/(sum(Q.values())), 'H2S': Q_S/(sum(Q.values())), }
# -------------------------------------------------------------------------------------------
plt.figure()
plt.plot(t_span, Q_M, label='CH4')
plt.plot(t_span, Q_C, label='CO2')
plt.plot(t_span, Q_S, label='Sulfur')
plt.legend()
plt.xlabel('Time [d]')
plt.ylabel('Gas Production Rate [m^3/d]')

plt.figure()
plt.plot(t_span, x_vol['CH4'], label='CH4')
plt.plot(t_span, x_vol['CO2'], label='CO2')
plt.plot(t_span, x_vol['H2S'], label='Sulfur')
plt.show()