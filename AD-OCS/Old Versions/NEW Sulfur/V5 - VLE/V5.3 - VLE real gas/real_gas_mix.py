# Solution for Anaerobic Digestion Thermodynamic Equilibrium - VLE
import numpy as np
from functions_GL import*
from scipy.optimize import fsolve
def f_VL_realgas(z,T,P):
    z = {
        'CH4': z[0],
        'CO2': z[1],
        'H2S': z[2],
        'H2O': z[3]
        } # mol/mol

    
    C = {# A, B, C, D, E, and F
        'CH4': [3.135E+1, -1.308E+3, 0, -3.261E+0, 2.942E-5, 2],
        'CO2': [1.336E+2, -4.735E+3, 0, -2.127E+1, 4.091E-2, 1],
        'H2S': [7.868E+1, -3.840E+3, 0, -1.120E+1, 1.885E-2, 1],
        'H2O': [6.593E+1, -7.228E+3, 0, -7.177E+0, 4.031E-6, 2]
    }

    Ant = {# A, B, C
        'CH4': [3.98950,   443.028,	 -0.49],
        'CO2': [6.81228,	1301.679,	-3.494],
        'H2S': [4.52887,	958.5870,	-0.539],
        'H2O': [6.20963,	2354.731,	 7.559]
    }


    kH = {# A, B, C and D
        'CH4': [-142.234, 26.5960, 32.308, -0.090],
        'CO2': [-183.691, 676.278, 39.068, -0.098],
        'H2S': [-186.884, 1090.12, 39.011, -0.094]
    }

    Tc = {
        'CH4': 190.4,
        'CO2': 304.1,
        'H2S': 373.2,
        'H2O': 647.3
    } # K

    Pc = {
        'CH4': 46.00*0.986923,
        'CO2': 73.80*0.986923,
        'H2S': 89.40*0.986923,
        'H2O': 221.2*0.986923
    } # atm

    w = {
        'CH4': 0.011,
        'CO2': 0.239,
        'H2S': 0.081,
        'H2O': 0.344
    } 

    # Evaluation
    Psat = {
        'CH4': Antoine(Ant['CH4'], T+273.15),
        'CO2': Antoine(Ant['CO2'], T+273.15),
        'H2S': Antoine(Ant['H2S'], T+273.15),
        'H2O': Antoine(Ant['H2O'], T+273.15)  
        } # atm

    H = {
        'CH4': Henry(kH['CH4'], T+273.15)*0.00986923,
        'CO2': Henry(kH['CO2'], T+273.15)*0.00986923,
        'H2S': Henry(kH['H2S'], T+273.15)*0.00986923  
        } # atm

    # Zed solution
    Zed = {
        'CH4': root_Zed(w['CH4'], Tc['CH4'], Pc['CH4'], P, T, H['CH4']),
        'CO2': root_Zed(w['CO2'], Tc['CO2'], Pc['CO2'], P, T, H['CO2']),
        'H2S': root_Zed(w['H2S'], Tc['H2S'], Pc['H2S'], P, T, H['H2S']),
        'H2O': root_Zed(w['H2O'], Tc['H2O'], Pc['H2O'], P, T, Psat['H2O']),
        }

    # Fugacity coefficient evaluation # Zed, w, Tc, Pc, P, T
    # a,b coefficients by mixing rules, args =  z, w, Tc, Pc, T

    meth_wat = mix_coeff(
        [z['CH4'], z['H2O']], 
        [w['CH4'], w['H2O']],
        [Tc['CH4'], Tc['H2O']],
        [Pc['CH4'], Pc['H2O']],
        T+273.15
        )

    carb_wat = mix_coeff(
        [z['CO2'], z['H2O']], 
        [w['CO2'], w['H2O']],
        [Tc['CO2'], Tc['H2O']],
        [Pc['CO2'], Pc['H2O']],
        T+273.15
        )

    sulf_wat = mix_coeff(
        [z['H2S'], z['H2O']], 
        [w['H2S'], w['H2O']],
        [Tc['H2S'], Tc['H2O']],
        [Pc['H2S'], Pc['H2O']],
        T+273.15
        )

    whole_mix = mix_coeff(
        [z['CH4'], z['CO2'], z['H2S'], z['H2O']], 
        [w['CH4'], w['CO2'], w['H2S'], w['H2O']],
        [Tc['CH4'], Tc['CO2'], Tc['H2S'], Tc['H2O']],
        [Pc['CH4'], Pc['CO2'], Pc['H2S'], Pc['H2O']],
        T+273.15    
    )

    mix_whole = {
        'a': whole_mix[0],
        'b': whole_mix[1]
    }

    amix = {
        'CH4': meth_wat[0],
        'CO2': carb_wat[0],
        'H2S': sulf_wat[0]
    }

    bmix = {
        'CH4': meth_wat[1],
        'CO2': carb_wat[1],
        'H2S': sulf_wat[1]
    }


    phi_V_mix = { # Zed, w, Tc, Pc, P, T, amix, bmix
        'CH4': fug_mix_PR(Zed['CH4'], w['CH4'], Tc['CH4'], Pc['CH4'], P, T+273.15, mix_whole['a'], mix_whole['b']),
        'CO2': fug_mix_PR(Zed['CO2'], w['CO2'], Tc['CO2'], Pc['CO2'], P, T+273.15, mix_whole['a'], mix_whole['b']),
        'H2S': fug_mix_PR(Zed['H2S'], w['H2S'], Tc['H2S'], Pc['H2S'], P, T+273.15, mix_whole['a'], mix_whole['b']),
        'H2O': 1
        }


    phi_V_sat = {
        'CH4': fug_coef_PR(Zed['CH4'], w['CH4'], Tc['CH4'], Pc['CH4'], Psat['CH4'], T+273.15),
        'CO2': fug_coef_PR(Zed['CO2'], w['CO2'], Tc['CO2'], Pc['CO2'], Psat['CO2'], T+273.15),
        'H2S': fug_coef_PR(Zed['H2S'], w['H2S'], Tc['H2S'], Pc['H2S'], Psat['H2S'], T+273.15),
        'H2O': 1
        }

    k = {
        'CH4': H['CH4']*phi_V_sat['CH4']/(P*phi_V_mix['CH4']),
        'CO2': H['CO2']*phi_V_sat['CO2']/(P*phi_V_mix['CO2']),
        'H2S': Psat['H2S']*phi_V_sat['H2S']/(P*phi_V_mix['H2S']),
        'H2O': Psat['H2O']*phi_V_sat['H2O']/(P*phi_V_mix['H2O'])    
    }

    # System Solution
    alpha0 = z['CH4'] + z['CO2'] + z['H2S']

    alpha = fsolve(
        func = RR, 
        x0   = alpha0, 
        args =(
            list(z.values()),
            list(k.values())
            )
        )
    alpha = alpha[0]
    k = [k['CH4'],k['CO2'], k['H2S'], k['H2O']]
    return alpha, k
#     # Visualization
# x = {
#     'CH4': z['CH4']/(1+alpha*(k['CH4']-1)),
#     'CO2': z['CO2']/(1+alpha*(k['CO2']-1)),
#     'H2S': z['H2S']/(1+alpha*(k['H2S']-1)),
#     'H2O': z['H2O']/(1+alpha*(k['H2O']-1)) 
#     } # mol/mol

# y = {
#     'CH4': x['CH4']*k['CH4'],
#     'CO2': x['CO2']*k['CO2'],
#     'H2S': x['H2S']*k['H2S'],
#     'H2O': x['H2O']*k['H2O'] 
#     } # mol/mol_wet

# V = alpha * F
# L = F - V

# Vy = {
#     'CH4': y['CH4']*V,
#     'CO2': y['CO2']*V,
#     'H2S': y['H2S']*V,
#     'H2O': y['H2O']*V 
#     } # kmol/d

# Lx = {
#     'CH4': x['CH4']*L,
#     'CO2': x['CO2']*L,
#     'H2S': x['H2S']*L,
#     'H2O': x['H2O']*L 
#     } # kmol/d


# y_dry = {
#     'CH4': y['CH4']/( sum(list(y.values())) - y['H2O'] ),
#     'CO2': y['CO2']/( sum(list(y.values())) - y['H2O'] ),
#     'H2S': y['H2S']/( sum(list(y.values())) - y['H2O'] )
#     } # mol/mol_dry 

# print('Mixture coefficients:')
# print('whole-mixture')
# pp(mix_whole)
# print('\nCH4-H2O')
# print('a', amix['CH4'])
# print('b', bmix['CH4'])
# print('\nCO2-H2O')
# print('a', amix['CO2'])
# print('b', bmix['CO2'])
# print('\nH2S-H2O')
# print('a', amix['H2S'])
# print('b', bmix['H2S'])
# print('')

# print('RESULTS')
# print('\nZed-factor:')
# pp(Zed)
# print('\nFugacity coeff mix - Vapour phase:')
# pp(phi_V_mix)
# print('\nFugacity coeff - Liquid phase:')
# pp(phi_V_sat)

# print('\nalpha:', alpha)
# print('alpha0:', alpha0)
# print('\nliquid molar fractions:')
# pp(x)
# print('\nvapour (wet basis) molar fractions:')
# pp(y)
# print('\nvapour (dry basis) molar fractions:')
# pp(y_dry)
# print('\nVapour flows [kmol/d]:')
# pp(Vy)
# print('\nLiquid flows [kmol/d]:')
# pp(Lx)
# print('\nFlows overlook: ')
# print('F = ', F, ' kmol/d')
# print('V = ', V, ' kmol/d')
# print('L = ', L, ' kmol/d')

    
