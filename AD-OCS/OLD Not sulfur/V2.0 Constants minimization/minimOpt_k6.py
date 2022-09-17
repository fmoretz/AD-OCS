# import scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, fsolve
from scipy.optimize import Bounds
from typing import List
from sklearn.metrics import r2_score
plt.close('all')

# import data and define the variables

from dataimport import *
from Identification import *
kd1 = kd[0]; kd2 = kd[1]
C_d = [0.1, 0.1]

D     = 1/HRT
mu_1 = alfa*D + C_d[0]*mu1_max

# create main function
def main():

    # define the bounds
    bd = Bounds([0,0],[20,20])

    # define the initial values
    gs = [200]

    # solve minimization problem
    y5exp  = q_M/X_2
    minsol = minimize(rmse, gs, args = (mu_2, y5exp), method = 'SLSQP', tol = 1e-10)
    y5model = y5eval(mu_2, minsol.x[0])

    # print out results
    fprint(minsol.x, minsol.success, minsol.nit, minsol.fun, [y5exp, y5model])

    # plot results
    x5 = mu_2
    regplot(x5, [y5exp, y5model], 'x = mu_2)', 'y = q_M/X_2', 'k6 regression', ['Experimental data', 'Model'], True)
    regplot(1/D, [y5exp, y5model], 'time [days]', 'y = q_M/X_2]', 'Model prediction vs Exp Data', ['Experimental data', 'Model'], True)


    k6x = np.linspace(0, 800, 20)
    rmsey = np.empty([len(k6x)])
    for i in range(0,len(k6x)):
            rmsey[i] = rmse([k6x[i]], mu_2, y5exp)

    plt.figure()
    plt.plot(k6x, rmsey)
    plt.title('RMSE for k6 regression')
    plt.xlabel('k6')
    plt.show()


    return 0;



# define the function for RSE
def rmse(x, mu_2, y5exp) -> float:
    y5list = []
    for i in range(0, len(XT)):
         y5list.append(y5eval(mu_2[i], x[0]))

    return np.sqrt(np.mean(y5exp-y5list)**2)

# define the function for the qC
def y5eval(mu_2, k6) -> float:
    return k6*mu_2

# utility function for printing and plotting
def fprint(sol, flag, itr, funcall, r2: List):
    print(f"\nSuccess: {flag}\n")
    print("==Results==")
    print(f"Solution: {sol}")
    print(f"k6: {sol}")
    print(f"RSE: {funcall}")
    print(f"Iterations: {itr}")
    print(f"R2 score: {r2_score(r2[0], r2[1])}")
    print("============")

def regplot(x, y: List, xlabel, ylabel, title, legend: List, showbool):
    plt.figure(figsize=(10,5))
    plt.plot(x, y[0], 'o', label = legend[0])
    plt.plot(x, y[1], '-', label = legend[1], linewidth=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.legend()
    plt.title(title)

    if showbool == True:
        plt.show()
    else:
        pass

if __name__ == "__main__":
    main()
