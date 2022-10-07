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
                print(t, t_change_loc, scale_loc)
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