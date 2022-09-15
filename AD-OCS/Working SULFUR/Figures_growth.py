'''Series of codes to try things related to sulfur modeling'''
import numpy as np
from main_V2 import*
import matplotlib.pyplot as plt
from PhysConstants import*

from plot import fitsubplotter 

# define gompertz function
def gompertz(x,a,b,c):
    return a*np.exp(-np.exp(b*np.exp(1)/a*(c-x) - 1))

def growth_SRB(x,a,b,c):                                                                    # Growth of SRB
    return b*np.exp((np.exp(1)*b/a*(c-x))-np.exp(np.exp(1)*b/a*(c-x)+1)+2)

'''1. Model represented in the main code'''
# Sulfur production model: dSs/dt = (1-y)/y * dXs/dt

Xs_max = np.zeros(len(t_span))
Xs = np.zeros(len(XT))
Ss = np.zeros(len(XT))
rho = np.zeros(len(XT))
rho_t = np.zeros(len(XT))
lam = np.zeros(len(XT))
mu_srb = np.zeros(len(XT))

for i in range(len(XT)):
    # Species differences
    Ss_max[i] = frac_sulfur*y_influent[i,4]*1000/64*S2[i]/y_influent[0,4]
    Xs_max[i] = Y_srb/(1-Y_srb)*Ss_max[i]   # maximum sulfur concentration g/L
    
    mu_srb[i] = (- X2[0] + X2[i])/(t_span[i] - t_span[0])     
    
    # Gompertz function
    Xs[i] = gompertz(t_span[i], Xs_max[i], mu_srb[i], 0)
    Ss[i] = Xs[i]*(1-Y_srb)/(Y_srb)
    rho[i] = growth_SRB(t_span[i], Xs_max[i], mu_srb[i], 0)


plt.figure(110)
plt.suptitle('As it is')
plt.subplot(3,1,1)
plt.plot(t_span, Xs, label="Gompertz")
# plt.plot(t_span, Xs_max, '--')
plt.ylabel('Conc [g/L]')
plt.grid(True)

plt.subplot(3,1,2)
plt.grid(True)
plt.plot(t_span, rho)

# plt.plot(t_span, mu_srb,'--')
plt.xlabel('Time [d]')
plt.ylabel('Growt [g/L/d]')

plt.subplot(3,1,3)
plt.grid(True)
plt.plot(t_span, Ss)
plt.xlabel('Time [d]')
plt.ylabel('Ss [mmol/L]')

'''2. Represent a Gompertz trend at each time step'''

size_iter = len(t_span)
Ss_max_ = np.zeros([size_iter])
Xs_max_ = np.zeros([size_iter])
Xs_ = np.zeros([size_iter,len(t_span)])
Ss_ = np.zeros([size_iter,len(t_span)])
rho_ = np.zeros([size_iter,len(t_span)])
rho_t = np.zeros([size_iter,len(t_span)])
lam = np.zeros([size_iter,len(t_span)])
mu_srb_ = np.zeros([size_iter,len(t_span)])

plt.figure(111)

for j in range(len(t_span)): 
    Ss_max_[j] = frac_sulfur*y_influent[j,4]*1000/64*S2[j]/y_influent[0,4]
    Xs_max_[j] = Y_srb/(1-Y_srb)*Ss_max[j]   # maximum sulfur concentration g/L
    
    for ii in range(len(t_span)):

        mu_srb_[j,ii] = (- X2[0] + X2[ii])/(t_span[ii] - t_span[0])     
    
        Xs_[j,ii]  = gompertz(t_span[ii], Xs_max_[j], mu_srb_[j,ii], 0)
        Ss_[j,ii]  = Xs_[j,ii]*(1-Y_srb)/(Y_srb)
        rho_[j,ii] = growth_SRB(t_span[ii], Xs_max_[j], mu_srb_[j,ii], 0)

    # get max of rho omitting nan
    print(np.nanmax(rho_))
    plt.suptitle('Gompertz trend at each time step')

    plt.subplot(2,1,1)
    # label with the index of the time step
 
    plt.plot(t_span, Xs_[j,:])
    # plt.plot(t_span, Xs_max, '--')
    plt.ylabel('Conc [g/L]')
    plt.grid(True)
    

    plt.subplot(2,1,2)
    plt.grid(True)
    plt.plot(t_span, rho_[j,:])
    plt.legend(ncol=2)
    # plt.plot(t_span, mu_srb,'--')
    plt.xlabel('Time [d]')
    plt.ylabel('Growth [g/L/d]')

plt.subplot(2,1,1)
plt.scatter(t_span, Xs, label="Gompertz")

'''3. Get the proper trend of growth rate - max of rho'''

size_iter = len(t_span)
Ss_max_ = np.zeros([size_iter])
Xs_max_ = np.zeros([size_iter])
Xs_ = np.zeros([size_iter,len(t_span)])
Ss_ = np.zeros([size_iter,len(t_span)])
rho_ = np.zeros([size_iter,len(t_span)])
rho_true = np.zeros([size_iter])
lam = np.zeros([size_iter,len(t_span)])
mu_srb_ = np.zeros([size_iter,len(t_span)])

plt.figure(112)
plt.suptitle('Gompertz trend at each time step and growth rate take as max')

for j in range(len(t_span)): 
    Ss_max_[j] = frac_sulfur*y_influent[j,4]*1000/64*S2[j]/y_influent[0,4]
    Xs_max_[j] = Y_srb/(1-Y_srb)*Ss_max[j]   # maximum sulfur concentration g/L
    
    for ii in range(len(t_span)):

        mu_srb_[j,ii] = (- X2[0] + X2[ii])/(t_span[ii] - t_span[0])     
    
        Xs_[j,ii]  = gompertz(t_span[ii], Xs_max_[j], mu_srb_[j,ii], 0)
        Ss_[j,ii]  = Xs_[j,ii]*(1-Y_srb)/(Y_srb)
        rho_[j,ii] = growth_SRB(t_span[ii], Xs_max_[j], mu_srb_[j,ii], 0)
    
    rho_true[j] = np.nanmax(rho_[j])
       # get max of rho omitting nan
      

    plt.subplot(2,1,1)
    # label with the index of the time step
 
    plt.plot(t_span, Xs_[j,:])
    # plt.plot(t_span, Xs_max, '--')
    plt.ylabel('Conc [g/L]')
    plt.grid(True)    

plt.subplot(2,1,2)
plt.grid(True)
plt.plot(t_span, rho_true)
plt.legend(ncol=2)
# plt.plot(t_span, mu_srb,'--')
plt.xlabel('Time [d]')
plt.ylabel('Growth [g/L/d]')

plt.subplot(2,1,1)
plt.scatter(t_span, Xs, color='red', label="Gompertz")
plt.plot(t_span, Xs, color = 'red', label="Gompertz")

plt.show()
