''' Contains the functions to be called in other codes'''

from dataimport import*

# Level dynamics function
def level_t(t, D, Qin, SR, h0, t_change):
    ''' 
    Function that returns the values of the liquid level at each time t.
    h0 [m] is the initial level
    t_change [h] is the time of the change in the influent flow - starts from zero
    D [1/d] - Dilution Rate
    Qin [m3/d] - Influent flow
    SR [m2] - Surface of the reactor
    '''
    t_loc = t-t_change # [h] - Time from the start of the change           
    return h0*np.exp(-D/24*t_loc) + (1-np.exp(-D/24*t_loc))*Qin/(D*SR)

def gompertz(x,a,b,c):                                                                      
    '''Gompertz function to model the SRB population'''
    return a*np.exp(-np.exp(b*np.exp(1)/a*(c-x) - 1))

def growth_SRB(x,a,b,c):                                                                    
    '''Function that returns the values of the SRB growth rate: it is the derivative of the Gompertz function'''
    return b*np.exp((np.exp(1)*b/a*(c-x))-np.exp(np.exp(1)*b/a*(c-x)+1)+2)

def f_deviations(t_span, t_change_vett, y_in_0):
    '''
    Function to be called in the main code.
    It returns the values of the deviated influents at each time t.
    Needs: 
    - t_span: time span of the simulation
    - t_change_vett: vector of the times at which the influent changes, defined by excel spreadsheet
    - y_in_0: vector of the initial values of the influent [S1in, S2in, Cin, Zin, XTin, Qin] and their measure units
    - Matrix T3: matrix of the deviations to be applied to the influent, derived from the Excel spreadsheet
    '''
    y_influent = np.zeros((len(t_span), len(y_in_0))) # Defines the y_in values which will be used in the main

    index = 0
    t_change_loc = t_change_vett[index]                # Get the time when the first deviation is applied
    scale_loc    = np.ones(len(y_in_0))    # Get the scale values for the first time change
   
    for i in range(len(t_span)):
        t = t_span[i]
        while True:
            if t < t_change_loc:        
                y_influent[i] = y_in_0*scale_loc
                print(round(t,2), t_change_loc, scale_loc)
                if i == (len(t_span)-1):  
                    break
                else:      
                    i += 1
                    t = t_span[i]
            else:        
                scale_loc    = T3.loc[t_change_loc].to_numpy()
                print(f'*** CHANGE at: {t_change_loc} ***')
                if index == (len(t_change_vett)-1):
                    t_change_loc = np.inf
                    print('Last Change done')
                else:
                    index = index+1        
                    t_change_loc = t_change_vett[index]
        return y_influent

def deviations_check(t, original, t_change_vett):
    '''Defines the deviation to be used in the ODE integration'''
    deviated = np.zeros(len(original)) 
    i = 0    
    while True: 
        if t < t_change_vett[0]: # if the time is before the first time change then the original values are used
            scale = 1
            break
        if t >= t_change_vett[-1]: # if the time is after the last time change then the last values are used
            t_change = t_change_vett[-1]
            scale = T3.loc[t_change].to_numpy()
            break
        if t > t_change_vett[i] and t < t_change_vett[i+1]:
            t_change = t_change_vett[i]
            scale = T3.loc[t_change].to_numpy()
            break
        else:
            i = i+1            
          
    deviated = original*scale           
    return deviated

def Antoine(x,T):
  '''
  Modified Antoine Equation
  lnPev = A + B/(T+C) + DlnT + E*T^F
  A, B, C, D, E, and F = fitted coefficients
  Pvap = saturation pressure in kPa
  T = temperature in K
  return: pressure of the species in kPa

  Normal Antoine Equation
  log10Pev = A - B/(T+C)
  A, B and C = fitted coefficients
  Pvap = saturation pressure in bar
  T = temperature in K
  return: pressure of the species in bar
  '''
  # return (np.exp( x[0] + x[1]/(T+x[2]) + x[3]*np.log(T) + x[4]*T**x[5] )) * 0.00986923
  return (10**( x[0] - x[1]/(T+x[2]) ) )* 0.986923

def Henry(x,T):
  '''
  Extended Henry's Equation
  lnH = A + B/T + ClnT + D*T
  A, B, C, and D = fitted coefficients
  H = henry's constants in kPa
  T = temperature in K
  return: henry's constants of the species in kPa
  '''
  return np.exp( x[0] + x[1]/T + x[2]*np.log(T) + x[3]*T )

def RR(alpha, z, k):
  '''
  Rachford-Rice equation, solution 
  for ideal gas/ideal liquid flash systems
  z: feed molar composition
  k: pressure ratio, defined as Psat/Ptot
  Psat: saturation pressure of the species
  Ptot: total pressure of the system
  return: after zero this function, it returns
  the value of alpha (vapor fraction)
  '''
  zdk = np.empty(shape=(len(z),))
  dk  = np.empty(shape=(len(z),))

  for i in range(len(z)):
    dk[i]  = (k[i]-1)
    zdk[i] = z[i]*dk[i]

  return  sum(zdk/(1+alpha*dk))

def root_Zed(w, Tc, Pc, P, T, Psat, display = False, u = 2, v = -1):
  '''
  Find the compressibility factor
  for a species with Peng-Robinson EOS.
  objective function: 
  0 = Zed^3 + alpha*Zed^2 + beta*Zed + gamma
  where:
  alpha = - 1 - B + u*B
  beta  = A + v*B^2 - u*B - u*B^2
  gamma = -A*B - v*B^2 - v*B^3
  S = 0.37464 + 1.5422 * w - 0.26992 * w**2
  k = ( 1 + S*(1-(T/Tc)**0.5) )**2
  a = (0.45724*k*(Rgas*Tc)**2)/Pc
  b = (0.07780*Rgas*Tc)/Pc
  A = (a*P)/(Rgas*T)**2
  B = (b*P)/(Rgas*T)
  return: value of Zed according to the pressure value 
  '''
  # Parameter evaluation
  Rgas = 0.082057338 # L*atm/K/mol
  S = 0.37464 + 1.54226 * w - 0.26992 * w**2
  k = ( 1 + S*(1-(T/Tc)**0.5) )**2
  a = (0.45724*k*(Rgas*Tc)**2)/Pc
  b = (0.07780*Rgas*Tc)/Pc
  A = (a*P)/((Rgas*T)**2)
  B = (b*P)/(Rgas*T)
  alpha = - 1 - B + u*B
  beta  = A + v*B**2 - u*B - u*B**2
  gamma = -A*B - v*B**2 - v*B**3

  # Solver coefficients for Peng-Robinson
  # Objective function: X^3 + p*X + q = 0 --> Zed = X - alpha/3
  p = beta - (alpha**2)/3
  q = 2*(alpha**3)/27 - alpha*beta/3 + gamma

  # Discrimant evaluation
  D = (q**2)/4 + (p**3)/27

  # Sign verification:
  if D > 0:
    X = ( -q/2 + (D**0.5) )**(1/3) + ( -q/2 - (D**0.5) )**(1/3)
    Zed_1 =  X - alpha/3
    Zed_2 = Zed_1
    Zed_3 = Zed_1
  elif D == 0:
    X = (q/2)**(1/3)
    Zed_1 = -2*X - alpha/3
    Zed_2 = X - alpha/3
    Zed_3 = Zed_2
  elif D < 0:
    r = (-p**3 / 27)**0.5
    theta = np.arccos( -q/(2*r) )
    Zed_1 = 2*r**(1/3)*np.cos(theta/3) - alpha/3
    Zed_2 = 2*r**(1/3)*np.cos( (2*m.pi + theta)/3 ) - alpha/3
    Zed_3 = 2*r**(1/3)*np.cos( (4*m.pi + theta)/3 ) - alpha/3
  
  Zed = [np.real(Zed_1), np.real(Zed_2), np.real(Zed_3)]
  
  if display == False:
    pass
  elif display == True:
    print('EOS roots:')
    print('Discriminant = ', D)
    print('Zed(1) = ', Zed[0])
    print('Zed(2) = ', Zed[1])
    print('Zed(3) = ', Zed[2])
    print('')
  
  # Choose the right Zed
  if P/Psat > 1:
    Zed_ = min(Zed) # Liquid phase
  else:
    Zed_ = max(Zed) # Gas phase

  return Zed_
  
def flash(u, z, Psat, H, P, phi_V_sat, phi_V):
  '''
  sysmple flash system 7x7 for 
  thermodynamic state solution
  '''

  # Unknown definition
  y = {
      'CH4': u[0],
      'CO2': u[1],
      'H2S': u[2],
      'H2O': 1-u[0]-u[1]-u[2] 
      } # mol/mol_wet
  
  x = {
      'CH4': u[3],
      'CO2': u[4],
      'H2S': u[5],
      'H2O': 1-u[3]-u[4]-u[5]
      } # mol/mol

  alpha = u[6]

  # Equation definition
  eq = np.empty(7)
  eq[0] = P*y['CH4']*phi_V['CH4'] - phi_V_sat['CH4']*Psat['CH4']*x['CH4']
  eq[1] = P*y['CO2']*phi_V['CO2'] - phi_V_sat['CO2']*Psat['CO2']*x['CO2']
  eq[2] = P*y['H2S']*phi_V['H2S'] - phi_V_sat['H2S']*Psat['H2S']*x['H2S']
  eq[3] = P*y['H2O']*phi_V['H2O'] - phi_V_sat['H2O']*Psat['H2O']*x['H2O']
  eq[4] = z['CH4'] - (1-alpha)*x['CH4'] - alpha*y['CH4']
  eq[5] = z['CO2'] - (1-alpha)*x['CO2'] - alpha*y['CO2']
  eq[6] = z['H2S'] - (1-alpha)*x['H2S'] - alpha*y['H2S']
  
  return eq

def fug_coef_PR(Zed, w, Tc, Pc, P, T):
    '''
    evaluation of the fugacity coeffcient
    from the Peng-Robinson EOS in a mixture.
    phi = exp(
      Bi/B*(Zed-1)-A/(2*2^0.5*B)*(Bi/B-2*(Ai/A)^0.5)*
      ln((Zed+B*(1+2^0.5))/(Zed+B*(1-2^0.5)))-ln(Zed-B))
    S = 0.37464 + 1.5422 * w - 0.26992 * w**2
    k = ( 1 + S*(1-(T/Tc)**0.5) )**2
    a = (0.45724*k*(Rgas*Tc)**2)/Pc
    b = (0.07780*Rgas*Tc)/Pc
    A = (a*P)/(Rgas*T)**2
    B = (b*P)/(Rgas*T)
    return: value of the fugacity coefficeint accordin to the system pressure
    '''
    # Parameter evaluation
    Rgas = 0.082057338 # L*atm/K/mol
    S = 0.37464 + 1.54226 * w - 0.26992 * w**2
    k = ( 1 + S*(1-(T/Tc)**0.5) )**2
    a = (0.45724*k*(Rgas*Tc)**2)/Pc
    b = (0.07780*Rgas*Tc)/Pc
    A = (a*P)/((Rgas*T)**2)
    B = (b*P)/(Rgas*T)

    # Parameter gathering
    c1 = Zed-1
    c2 = A/(2*(2**0.5)*B)
    c3 = np.log( (Zed+B*(1+2**0.5))/(Zed+B*(1-2**0.5)) )
    c4 = np.log(Zed-B)
    phi = np.exp( c1 - c2 * c3 - c4 )

    # check the value
    if np.isnan(phi) == True:
      phi = 0
    else:
      pass
    return phi
  
def mix_coeff(z, w, Tc, Pc, T):
  
  # Parameter evaluation
  Rgas = 0.082057338 # L*atm/K/mol

  ai  = np.empty(shape=(len(z)))
  bi  = np.empty(shape=(len(z)))
  S   = np.empty(shape=(len(z)))
  k   = np.empty(shape=(len(z)))
  zai = np.empty(shape=(len(z)))

  for i in range(len(z)):
    S[i]   = 0.37464 + 1.54226 * w[i] - 0.26992 * w[i]**2
    k[i]   = ( 1 + S[i]*(1-(T/Tc[i])**0.5) )**2
    ai[i]  = (0.45724*k[i]*(Rgas*Tc[i])**2)/Pc[i]
    bi[i]  = (0.07780*Rgas*Tc[i])/Pc[i]
    zai[i] = z[i]*(ai[i]**0.5)

  a = sum(zai)**2
  b = sum(bi)/len(bi)

  return [a, b]

def fug_mix_PR(Zed, w, Tc, Pc, P, T, amix, bmix):
  '''
  evaluation of the fugacity mixture coeffcient
  from the Peng-Robinson EOS in a mixture.
  phi = exp(
    Bi/B*(Zed-1)-A/(2*2^0.5*B)*(Bi/B-2*(Ai/A)^0.5)*
    ln((Zed+B*(1+2^0.5))/(Zed+B*(1-2^0.5)))-ln(Zed-B))
  S = 0.37464 + 1.5422 * w - 0.26992 * w**2
  k = ( 1 + S*(1-(T/Tc)**0.5) )**2
  ai = (0.45724*k*(Rgas*Tc)**2)/Pc
  bi = (0.07780*Rgas*Tc)/Pc
  aij = (ai*aj)*(1-kij) = (sum(zi*(ai)**0.5))**2
  bij = (bi + bj)/2
  A = (a*P)/(Rgas*T)**2
  B = (b*P)/(Rgas*T)
  return: value of the fugacity coefficeint accordin to the system pressure
  '''
  # Parameter evaluation
  Rgas = 0.082057338 # L*atm/K/mol
  S = 0.37464 + 1.54226 * w - 0.26992 * w**2
  k = ( 1 + S*(1-(T/Tc)**0.5) )**2
  a = (0.45724*k*(Rgas*Tc)**2)/Pc
  b = (0.07780*Rgas*Tc)/Pc
  Amix = (amix*P)/((Rgas*T)**2)
  Bmix = (bmix*P)/(Rgas*T)
  A = (a*P)/((Rgas*T)**2)
  B = (b*P)/(Rgas*T)

  # Parameter gathering
  c1 = B/Bmix*(Zed-1)
  c2 = A/(2*(2**0.5)*B)
  c3 = B/Bmix - 2*(A/Amix)**0.5
  c4 = np.log( (Zed+Bmix*(1+2**0.5))/(Zed+Bmix*(1-2**0.5)) )
  c5 = np.log(Zed-Bmix)
  phi_mix = np.exp( c1 + c2 * c3 * c4 - c5)

  # check the value
  if np.isnan(phi_mix) == True:
    phi = 0
  else:
    pass
  return phi_mix


def headspace_dynamics_discr(Nin, Nin_O2, P, T, Vgas, t):
  # Nin = [CH4, CO2, H2S, H2O]
    KI_O2 = 1.9e-7 # [kg/m3]
    dt = (t[-1] - t[0])/len(t)
    k    = 0.6
    alfa = 1
    beta = 0.1
    
    R = 0.0821 # L * atm/mol/K
    

    n_in = {
        'H2S': Nin[:,2],
        'O2':  Nin_O2,
        'SX':  0.000,
        'H2O': Nin[:,3]
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

    n_out_inhib = np.zeros(len(t))
    I_O2 = np.zeros(len(t))
    

    n_out['H2S'][0] = n_in['H2S'][0]
    n_out['O2'][0]  = n_in['O2'][0]
    n_out['H2O'][0] = n_in['H2O'][0]

    x_out['H2S'][0] = n_in['H2S'][0] / (n_out['H2O'][0] + n_out['H2S'][0] + n_out['O2'][0] + Nin[0,0] + Nin[0,1])
    x_out['O2'][0]  = n_in['O2'][0]  / (n_out['H2O'][0] + n_out['H2S'][0] + n_out['O2'][0] + Nin[0,0] + Nin[0,1])

    c_out['H2S'][0] = x_out['H2S'][0] * P/R/(T+273.15) # mol/L
    c_out['O2'][0]  = x_out['O2'][0]  * P/R/(T+273.15) # mol/L

    m_out['H2S'][0] = c_out['H2S'][0] * PM['H2S'] * 1000 # g/m3
    m_out['O2'][0]  = c_out['O2'][0]  * PM['O2'] * 1000 # g/m3

    I_O2[0] = 1/(1+m_out['O2'][0]/KI_O2/1000)
    n_out_inhib[0] = Nin[0,0]*I_O2[0]
    r_SOB = np.zeros(len(t))
    r_SOB[0] = (k * m_out['H2S'][0]**alfa * m_out['O2'][0]**beta) # gS/m3/h
     

    for i in range(1, len(t)):
                             #                                                           CH4       CO2
        x_out['H2S'][i-1] = n_out['H2S'][i-1] / (n_out['H2O'][i-1] + n_out['H2S'][i-1] + n_out['O2'][i-1] + Nin[i,0] + Nin[i,1])
        x_out['O2'][i-1]  = n_out['O2'][i-1]  / (n_out['H2O'][i-1] + n_out['H2S'][i-1] + n_out['O2'][i-1] + Nin[i,0] + Nin[i,1])
        
        c_out['H2S'][i-1] = x_out['H2S'][i-1] * P/R/(T+273.15) # mol/L
        c_out['O2'][i-1]  = x_out['O2'][i-1]  * P/R/(T+273.15) # mol/L
        
        m_out['H2S'][i-1] = c_out['H2S'][i-1] * PM['H2S'] * 1000 # g/m3
        m_out['O2'][i-1]  = c_out['O2'][i-1]  * PM['O2'] * 1000 # g/m3

        r_SOB[i] = (k * m_out['H2S'][i-1]**alfa * m_out['O2'][i-1]**beta) # gS/m3/h

        I_O2[i] = 1/(1+m_out['O2'][i-1]/KI_O2/1000) 
        n_out_inhib[i] = Nin[i,0]*I_O2[i]

        n_out['H2S'][i] = (n_in['H2S'][i]*dt + (n_out['H2S'][i-1])  - r_SOB[i]*dt* Vgas[i] / PM['H2S'])/(1+dt)
        n_out['O2'][i]  = (n_in['O2'][i]*dt  + (n_out['O2'][i-1])   - r_SOB[i]*dt* Vgas[i] / (2 * PM['O2']))/(1+dt)
        n_out['SX'][i]  = (n_in['SX']*dt  + (n_out['SX'][i-1])   + r_SOB[i]*dt* Vgas[i] / PM['SX'])/(1+dt)
        
        n_out['H2O'][i] = n_in['H2O'][i] + (n_out['H2S'][i-1] - n_in['H2S'][i]) # n_in['H2O'] + nu['H2O'] * (k * m_out['H2S'][i]**alfa * m_out['O2'][i]**beta)*( dt/(1+dt) ) * Vgas /
    
    return n_out, r_SOB, n_out_inhib, I_O2



def Headspace_reactions(y, RC, k_SOB, alfa, beta, V_gas, N_V):
        
        n_S =  y[0]                             # [mol/h] - Sulfur flowrate
        n_O2 = y[1]                              # [mol/h] - Oxygen flowrate  
      
        n_out = N_V[0]+N_V[1]+N_V[i,3] + n_S + n_O2 # [mol/h] - Total Outflow
        
        y_S_in = n_S_in/(sum(N_V[:])+n_O2_in)         # [-] - Sulfur molar fraction in the inflow
        y_O_in = n_O2_in/(sum(N_V[:])+n_O2_in)        # [-] - Sulfur molar fraction in the inflow

        y_S_loc = n_S/n_out                     # [-] - Sulfur fraction in the headspace
        y_O2_loc = n_O2/n_out                   # [-] - Oxygen fraction in the headspace       
        
        C_S  = y_S_loc*P_dig/(Rgas_L_atm_K*T)   # [mol/L] - Sulfur concentration in gas phase atm/(L*atm/mol/K*L) = mol/L
        C_O2 = y_O2_loc*P_dig/(Rgas_L_atm_K*T)  # [mol/L] - Oxygen concentration in gas phase

        S_S = C_S*MW_S*1e+3                     # [g/m3] - Sulfur concentration in gas phase
        S_O2 = C_O2*MW_O*1e+3                   # [g/m3] - Oxygen concentration in gas phase

        r_SOB = k_SOB*(S_S**alfa)*(S_O2**beta)                     # [gS/m3/h] - Reaction rate SOB
       
        
        
        BM_S  = n_S_in     - n_S  - V_gas*r_SOB/MW_S            # [mol/h]   - Sulfur Balance
        BM_O2 = n_O2_in    - n_O2 - V_gas*r_SOB/RC/MW_O         # [mol/h]   - Oxygen Balance
        
        if BM_S and BM_O2 < 1e-10:
            r_vett[i] = r_SOB/60                                    # [mg/L/min] - Reaction rate
            S_vett[0] = S_S                                       # [mg/L] - Sulfur concentration
            S_vett[1] = y_S_in*P_dig/(Rgas_L_atm_K*T)*MW_S*1e+3   # [mg/L] - Influent conc of Sulfur
            S_vett[2] = S_O2                                      # [mg/L] - O2 concentration
            S_vett[3] = y_O_in*P_dig/(Rgas_L_atm_K*T)*MW_O*1e+3   # [mg/L] - Influent conc of Sulfur
            print('r = ', r_SOB/60,'[mg/L/min]')
            print('S_S,in = ', y_S_in*P_dig/(Rgas_L_atm_K*T)*MW_S*1e+3,'[mg/L]', 'S_O2,in', y_O_in*P_dig/(Rgas_L_atm_K*T)*MW_O*1e+3,'[mg/L]')
            print('S_S = ', S_S,'[mg/L]', 'S_O2 = ', S_O2,'[mg/L]')
        return [BM_S, BM_O2]