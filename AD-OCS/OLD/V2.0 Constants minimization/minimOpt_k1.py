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
    gs = [20]

    # solve minimization problem
    y4exp = D*(S1_in - S1) + k_hyd*XT
    minsol = minimize(rmse, gs, args = (X_1, mu_1, y4exp), method = 'SLSQP', tol = 1e-10)
    y4model = y4eval(X_1, mu_1, minsol.x[0])

    # print out results
    fprint(minsol.x, minsol.success, minsol.nit, minsol.fun, [y4exp, y4model])

    # plot results
    x4 = X_1*mu_1
    regplot(x4, [y4exp, y4model], 'x = X_1*(alpha*D + kd1)', 'y = D*(S1in-S1) + k_hyd*XT , 'k1 regression', ['Experimental data', 'Model'], True)
    regplot(1/D, [y4exp, y4model], 'time [days]', 'y = D*(S1in-S1) + k_hyd*XT [mmol/L/d]', 'Model prediction vs Exp Data', ['Experimental data', 'Model'], True)


    k1x = np.linspace(0, 40, 40)
    rmsey = np.empty([len(k1x)])
    for i in range(0,len(k1x)):
            rmsey[i] = rmse([k1x[i]], X_1, mu_1, y4exp)

    plt.figure()
    plt.plot(k1x, rmsey)
    plt.title('RMSE for k1 regression')
    plt.xlabel('k1')
    plt.show()


    return 0;



# define the function for RSE
def rmse(x, X_1, mu_1, y4exp) -> float:
    y4list = []
    for i in range(0, len(D)):
         y4list.append(y4eval(X_1[i], mu_1[i], x[0]))

    return np.sqrt(np.mean(y4exp-y4list)**2)

# define the function for the qC
def y4eval(X_1, mu_1, k1) -> float:
    return k1*X_1*mu_1

# utility function for printing and plotting
def fprint(sol, flag, itr, funcall, r2: List):
    print(f"\nSuccess: {flag}\n")
    print("==Results==")
    print(f"Solution: {sol}")
    print(f"k1: {sol}")
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

def surfplot(x, y, z, xlabel, ylabel, zlabel, title, showbool):
    import mpl_toolkits.mplot3d.axes3d as axes3d
    from matplotlib import cm
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X, Y, z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_title(title)

    if showbool == True:
        plt.show()
    else:
        pass



if __name__ == "__main__":
    main()
