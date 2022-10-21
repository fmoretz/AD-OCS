import numpy as np
from matplotlib import pyplot as plt
import time
from pprint import pprint as pp

def headspace_dynamics_discr(Nin, P, T, Vgas, tau_head):

    t = np.linspace(0, tau_head, 1000)
    k    = 0.6
    alfa = 1
    beta = 0.1
    
    R = 0.0821 # L * atm/mol/K
    

    n_in = {
        'H2S': Nin[0],
        'O2':  Nin[2],
        'SX':  0.000,
        'H2O': Nin[1]
    } # mol/hours

    PM = {
        'H2S': 34.1,
        'O2':  36.0,
        'SX':  32.1,
        'H2O': 18.0
    }
    

    n_out = {
        'H2S': np.zeros(len(t)),
        'O2':  np.zeros(len(t)),
        'SX':  np.zeros(len(t)),
        'H2O': np.zeros(len(t))
    } # mol/hours

    c_out = {
        'H2S': np.zeros(len(t)),
        'O2':  np.zeros(len(t))
    } # mol/hours

    m_out = {
        'H2S': np.zeros(len(t)),
        'O2':  np.zeros(len(t))
    } # mol/hours

    x_out = {
        'H2S': np.zeros(len(t)),
        'O2':  np.zeros(len(t))
    } # mol/hours

    n_out['H2S'][0] = n_in['H2S']
    n_out['O2'][0]  = n_in['O2']
    n_out['H2O'][0] = n_in['H2O']

    x_out['H2S'][0] = n_in['H2S'] / (sum(n_in.values()) - n_in['SX'] + 2537.39 + 706.49)
    x_out['O2'][0]  = n_in['O2']  / (sum(n_in.values()) - n_in['SX'] + 2537.39 + 706.49)

    c_out['H2S'][0] = x_out['H2S'][0] * P/R/T # mol/L
    c_out['O2'][0]  = x_out['O2'][0]  * P/R/T # mol/L

    m_out['H2S'][0] = c_out['H2S'][0] * PM['H2S'] * 1000 # g/m3
    m_out['O2'][0]  = c_out['O2'][0]  * PM['O2'] * 1000 # g/m3

    print('Simulation is on going...')
    t0 = time.time()

    for i in range(1, len(t)):
    #                                                                                     CH4       CO2
        x_out['H2S'][i-1] = n_out['H2S'][i-1] / (n_out['H2O'][i-1] + n_out['H2S'][i-1] + n_out['O2'][i-1] + 2537.39 + 706.49)
        x_out['O2'][i-1]  = n_out['O2'][i-1]  / (n_out['H2O'][i-1] + n_out['H2S'][i-1] + n_out['O2'][i-1] + 2537.39 + 706.49)
        
        c_out['H2S'][i-1] = x_out['H2S'][i-1] * P/R/T # mol/L
        c_out['O2'][i-1]  = x_out['O2'][i-1]  * P/R/T # mol/L
        
        m_out['H2S'][i-1] = c_out['H2S'][i-1] * PM['H2S'] * 1000 # g/m3
        m_out['O2'][i-1]  = c_out['O2'][i-1]  * PM['O2'] * 1000 # g/m3
        
        n_out['H2S'][i] = (n_in['H2S']*dt + (n_out['H2S'][i-1])  - (k * m_out['H2S'][i-1]**alfa * m_out['O2'][i-1]**beta)*dt* Vgas / PM['H2S'])/(1+dt)
        n_out['O2'][i]  = (n_in['O2']*dt  + (n_out['O2'][i-1])  - (k * m_out['H2S'][i-1]**alfa * m_out['O2'][i-1]**beta)*dt* Vgas / (2 * PM['O2']))/(1+dt)
        n_out['SX'][i]  = (n_in['SX']*dt  + (n_out['SX'][i-1])  + (k * m_out['H2S'][i-1]**alfa * m_out['O2'][i-1]**beta)* dt* Vgas / PM['SX'])/(1+dt)
        
        n_out['H2O'][i] = n_in['H2O'] + (n_out['H2S'][i-1] - n_in['H2S']) # n_in['H2O'] + nu['H2O'] * (k * m_out['H2S'][i]**alfa * m_out['O2'][i]**beta)*( dt/(1+dt) ) * Vgas /

    
    return [n_out['H2S'], n_out['O2'], n_out['SX'], n_out['H2O']]



