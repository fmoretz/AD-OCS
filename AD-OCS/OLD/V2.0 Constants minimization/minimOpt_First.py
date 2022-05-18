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
S1_in  = float(T2["S1in"])
S2_in  = float(T2["S2in"])
D     = 1/T1["HRT"]
S_1   = T1["S1"]
S_2   = T1["S2"]
XT    = T1["XT"]
q_M   = T1["q_CH4"]        # Experimental points for the regression and comparison

# create main function
def main():

    # define the bounds
    bd = Bounds([0,0],[1000,1000])

    # define the initial values
    gs = [1,50]

    fargs = (D, S1, S2, XT, S1_in, S2_in, k_hyd)

    # solve minimization problem
    minsol = minimize(rmse, gs, args = fargs, method = 'SLSQP', bounds = bd, tol = 1e-10)
    qMmodel = qMeval(D, S1, S2, XT, S1_in, S2_in, k_hyd, minsol.x[0], minsol.x[1]) # output are AA and BB

    # print out results
    fprint(minsol.x, minsol.success, minsol.nit, minsol.fun, [q_M, qMmodel])

    # plot results
    regplot(1/D, [q_M, qMmodel], 'time [days]', 'qM [mmol/L/d]', 'Methane production regression', ['Experimental data', 'Model'], False)

    # plot function surface
    AAx = np.linspace(0, 4, 40)
    BBx = np.linspace(0, 10, 40)


    rmsez = np.empty([len(AAx), len(BBx)])
    for i in range(0,len(AAx)):
        for j in range(0,len(BBx)):
            rmsez[i][j] = rmse([AAx[i], BBx[j]], D, S1, S2, XT, S1_in, S2_in, k_hyd)

    surfplot(AAx, BBx, rmsez, 'AA', 'BB', 'rmse', 'rmse over AA and BB values', True)

    k2x = np.linspace(0.01, 2000, 20)
    k3x = np.linspace(50, 2000, 20)

    rmsez = np.empty([len(k2x), len(k3x)])
    for i in range(0,len(k2x)):
        for j in range(0,len(k3x)):
            rmsez[i][j] = rmse([k6/k3x[i], k2x[j]/k1], D, S1, S2, XT, S1_in, S2_in, k_hyd)

    surfplot(k2x, k3x, rmsez, 'k2', 'k3', 'rmse', 'rmse over k2 and k3 values', True)

    return 0;



# define the function for RSE - minimization on AA and BB
def rmse(x,D, S1, S2, XT, S1_in, S2_in, k_hyd) -> float:
    qMlist = []
    for i in range(0, len(D)):
         qMlist.append(qMeval(D[i], S1[i], S2[i], XT[i], S1_in, S2_in, k_hyd, x[0], x[1]))

    return np.sqrt(np.mean(q_M-qMlist)**2)

# define the function for the qM
def qMeval(D, S1, S2, XT, S1_in, S2_in, k_hyd, AA, BB):
    return AA*D*(S2_in - S2) + AA*BB*((S1_in-S1) + k_hyd*XT)

# utility function for printing and plotting
def fprint(sol, flag, itr, funcall, r2: List):
    print(f"\nSuccess: {flag}\n")
    print("==Results==")
    print(f"Solution: {sol}")
    print(f"k2: {sol[1]*k1}")
    print(f"k3: {k6/sol[0]}")
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
