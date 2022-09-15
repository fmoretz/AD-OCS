''' Functions to be called in the main code '''

from Influent import*

def gompertz(x,a,b,c):                                                                      # Gompertz function
    return a*np.exp(-np.exp(b*np.exp(1)/a*(c-x) - 1))

def growth_SRB(x,a,b,c):                                                                    # Growth of SRB
    return b*np.exp((np.exp(1)*b/a*(c-x))-np.exp(np.exp(1)*b/a*(c-x)+1)+2)

def AD_OCS_Model(x,t,alfa,mu_max,Ks,KI2,KH,Pt,kLa,D,k,kd,N_bac,N_S1,X2_0,t_0,y_in,t_change_vett):   # Model deviations function
    '''Main function used in the code: represents the AD-OCS model equations to be integrated'''   
    XT, X1, X2, Z, S1, S2, S2_new, C = x
    
    
    XT_in_0 = y_in_0[4] # Initial value of the influent particulate

    y_influent = deviations_check(t, y_in, t_change_vett)
    
    S1in = y_influent[0]      # [gCOD/L]
    S2in = y_influent[1]      # [mmol/L]
    Cin  = y_influent[2]      # [mmol/L]
    Zin  = y_influent[3]      # [mmol/L]
    XTin = y_influent[4]      # [gCOD/L]
   
    mu1 = mu_max[0]*(S1/(S1+Ks[0]))                                                                  # Monod
    mu2 = mu_max[1]*(S2/(S2+Ks[1]+S2**2/KI2))                                                        # Haldane

    qM  = k[5]*mu2*X2                                                                                # [mmol/L/day] - Methane molar flow 
    CO2 = C + S2 - Z                                                                                 # [mmol/L]     - CO2 Dissolved
    phi = CO2 + KH*Pt + qM/kLa
    Pc  = (phi - (phi**2- 4*KH*Pt*CO2)**0.5)/(2*KH)                                                  # [atm] - Partial pressure CO2
    qC  = kLa*(CO2 - KH*Pc)                                                                          # [mmol/L/d] - Carbon Molar Flow

    Ss_max = 0.02*XT_in*1000/64*S2/XT_in_0
    Xs_max = Y_srb/(1-Y_srb)*Ss_max   # maximum sulfur concentration g/L
    if t > 0:
        mu_max_srb = np.nan_to_num((X2 - X2_0)/(t - t_0), nan=0, neginf=0)        
    else:
        mu_max_srb = 0
    rho_srb = growth_SRB(t, Xs_max, mu_max_srb, 0)  
        
    dXT = D*(XTin - XT) - k[6]*XT                                                                    # Evolution of particulate
    dX1 = (mu1 - alfa*D - kd[0])*X1                                                                  # Evolution of biomass 1 (acidogen.)
    dX2 = (mu2 - alfa*D - kd[1])*X2                                                                  # Evolution of biomass 2 (methanogen)
    dZ  = D*(Zin - Z) + (k[0]*N_S1 - N_bac)*mu1*X1 - N_bac*mu2*X2 + kd[0]*N_bac*X1 + kd[1]*N_bac*X2  # Evolution of alcalinity;
    dS1 = D*(S1in - S1) - k[0]*mu1*X1 + k[6]*XT                                                      # Evolution of organic substrate
    dS2 = D*(S2in - S2) + k[1]*mu1*X1 - k[2]*mu2*X2                                                  # Evolution of VFA
    dS2_new = D*(S2in - S2_new) + k[1]*mu1*X1 - k[2]*mu2*X2 - rho_srb/Y_srb                              # Evolution of VFA - affected by sulfur
    dC  = D*(Cin - C)   + k[3]*mu1*X1 + k[4]*mu2*X2 - qC                                             # Evolution of inorganic carbon

    dxdt = [dXT, dX1, dX2, dZ, dS1, dS2, dS2_new,  dC]

    return dxdt

def f_deviations(t_span, t_change_vett, y_in_0):
    '''
    Function to be called in the main code.
    It returns the values of the deviated influents at each time t.
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