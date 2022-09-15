import numpy as np
from main import*
from matplotlib import pyplot as plt

# define gompertz function
def gompertz(x,a,b,c):
    return a*np.exp(-np.exp(b*np.exp(1)/a*(c-x) - 1))

def growth_SRB(x,a,b,c):                                                                    # Growth of SRB
    return b*np.exp((np.exp(1)*b/a*(c-x))-np.exp(np.exp(1)*b/a*(c-x)+1)+2)



Xs = np.zeros(len(XT))
Ss = np.zeros(len(XT))
rho = np.zeros(len(XT))
lam = np.zeros(len(XT))
mu_srb = np.zeros(len(XT))

print(Xs_max)
print(mu_srb)

plt.figure(105)


for i in range(len(XT)):
    Xs_max = 0.047  # maximum sulfur concentration g/L
    mu_srb[i] = (- X2[0] + X2[i])/(t_span[i] - t_span[0]) 

    Xs[i] = gompertz(t_span[i], Xs_max, mu_srb[i],0)
    rho[i] = growth_SRB(t_span[i], Xs_max, mu_srb[i],0)

plt.plot(t_span, Xs, label="Xs")
plt.grid(True)
plt.show()

# Xs = np.zeros([len(Xs_max),len(t_span)])
# for j in range(len(Xs_max)):
#     Xs_loc = Xs_max[j]
#     mu_srb_loc = mu_srb[j]
#     print(j)
#     for i in range(len(t_span)): # evaluate gompertz function
#         Xs[j,i] = gompertz(t_span[i], Xs_max[j], mu_srb[j], 0)

#     plt.plot(t_span, Xs[j], label=j)
# plt.legend(ncol=4)


