'''Here are present the model equations to be integrated:
    AMOCO_HN: AM2 with modifications accounting for hydrolysis and nitrgoen kinetics: Ficara et al. (2012)
    AD_OCS: AMOCOHN with kinetics for SRB and sulfur modelling. Own work.'''
from dataimport import*
from PhysConstants  import*
from functions import deviations_check, gompertz

def AD_OCS(x,t,alfa,mu_max,Ks,KI2,KH,Pt,kLa,D,k,kd,N_bac,N_S1,X2_0,t_0,y_in,t_change_vett, t_span):   # Model deviations function
    '''Main function used in the code: represents the AD-OCS model equations to be integrated'''   
    XT, X1, X2, Z, S1, S2, C = x
    
    XT_in_0 = y_in_0[4]         # Initial value of the influent particulate
    
    y_influent = deviations_check(t, y_in, t_change_vett) # Evaluate the influent deviations with an external function
    
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

    Ss_max = 0.02*XT_in*1000/64*S2/XT_in_0      # maximum sulfur concentration [g/L] - 0.02 is the maximum sulfur concentration in the influent [gCOD/gCOD]
    Xs_max = Y_srb/(1-Y_srb)*Ss_max             # maximum biomass concentration [gCOD/L] - by realtionship between Ss and Xs (microbioloigy)
    
    mu_srb = (+ X2 - X2_0)/max(1e-14,t-t_0)      # [g/L/d]    - Gompertz parameter for SRB growth
    Xs  = gompertz(t, Xs_max, mu_srb, 0)           # [g/L]      - Sulfate Reducing Bacteria - Gompertz
    rho_srb = mu2*Xs/X2 

        
    print('t:',round(t,4))
    dXT = D*(XTin - XT) - k[6]*XT                                                                    # Evolution of particulate
    dX1 = (mu1 - alfa*D - kd[0])*X1                                                                  # Evolution of biomass 1 (acidogen.)
    dX2 = (mu2 - alfa*D - kd[1])*X2                                                                  # Evolution of biomass 2 (methanogen)
    dZ  = D*(Zin - Z) + (k[0]*N_S1 - N_bac)*mu1*X1 - N_bac*mu2*X2 + kd[0]*N_bac*X1 + kd[1]*N_bac*X2  # Evolution of alcalinity;
    dS1 = D*(S1in - S1) - k[0]*mu1*X1 + k[6]*XT                                                      # Evolution of organic substrate
    dS2 = D*(S2in - S2) + k[1]*mu1*X1 - k[2]*mu2*X2 - rho_srb/Y_srb                                  # Evolution of VFA
    dC  = D*(Cin - C)   + k[3]*mu1*X1 + k[4]*mu2*X2 - qC                                             # Evolution of inorganic carbon
    
    dxdt = [dXT, dX1, dX2, dZ, dS1, dS2, dC]
    
    return dxdt

def AMOCO_HN(x,t,alfa,mu_max,Ks,KI2,KH,Pt,kLa,D,k,kd,N_bac,N_S1,X2_0,t_0,y_in,t_change_vett):   # Model deviations function
    '''Represents the AMOCO_HN model equations to be integrated'''   
    XT, X1, X2, Z, S1, S2, C = x
   
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
          
    dXT = D*(XTin - XT) - k[6]*XT                                                                    # Evolution of particulate
    dX1 = (mu1 - alfa*D - kd[0])*X1                                                                  # Evolution of biomass 1 (acidogen.)
    dX2 = (mu2 - alfa*D - kd[1])*X2                                                                  # Evolution of biomass 2 (methanogen)
    dZ  = D*(Zin - Z) + (k[0]*N_S1 - N_bac)*mu1*X1 - N_bac*mu2*X2 + kd[0]*N_bac*X1 + kd[1]*N_bac*X2  # Evolution of alcalinity;
    dS1 = D*(S1in - S1) - k[0]*mu1*X1 + k[6]*XT                                                      # Evolution of organic substrate
    dS2 = D*(S2in - S2) + k[1]*mu1*X1 - k[2]*mu2*X2                                                  # Evolution of VFA
    dC  = D*(Cin - C)   + k[3]*mu1*X1 + k[4]*mu2*X2 - qC                                             # Evolution of inorganic carbon

    dxdt = [dXT, dX1, dX2, dZ, dS1, dS2, dC]
    return dxdt
