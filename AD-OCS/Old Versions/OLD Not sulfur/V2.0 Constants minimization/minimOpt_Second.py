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
C_in  = float(T2["Cin"])
S1_in = float(T2["S1in"])
D     = 1/T1["HRT"]
S1   = T1["S1"]
XT   = T1["XT"]
C    = T1["C"]
q_C  = T1["q_C"]
q_M  = T1["q_CH4"]         # Experimental points for the regression and comparison

# create main function
def main():

    # define the bounds
    bd = Bounds([0,0],[20,20])

    # define the initial values
    gs = [1,50]

    # solve minimization problem
    minsol = minimize(rmse, gs, args = (D, S1, XT, q_M, C, C_in, S1_in, k_hyd), method = 'SLSQP', bounds = bd, tol = 1e-10)
    qCmodel = qCeval(D, S1, XT, q_M, C, C_in, S1_in, k_hyd, minsol.x[0], minsol.x[1])

    # print out results
    fprint(minsol.x, minsol.success, minsol.nit, minsol.fun, [q_C, qCmodel])

    # plot results
    regplot(1/D, [q_C, qCmodel], 'time [days]', 'q_C [mmol/L/d]', 'Carbon dioxide production regression', ['Experimental data', 'Model'], False)

    # plot function surface
    CCx = np.linspace(0, 4, 40)
    DDx = np.linspace(0, 10, 40)


    rmsez = np.empty([len(CCx), len(DDx)])
    for i in range(0,len(CCx)):
        for j in range(0,len(DDx)):
            rmsez[i][j] = rmse([CCx[i], DDx[j]], D, S1, XT, q_M, C, C_in,  S1_in, k_hyd)

    surfplot(CCx, DDx, rmsez, 'CC', 'DD', 'rmse', 'rmse over CC and DD values', False)

    k4x = np.linspace(0, 200, 20)
    k5x = np.linspace(0, 200, 20)

    rmsez = np.empty([len(k4x), len(k5x)])
    for i in range(0,len(k4x)):
        for j in range(0,len(k5x)):
            rmsez[i][j] = rmse([k4x[i]/k1, k5x[j]/k6], D, S1, XT, q_M, C, C_in, S1_in, k_hyd)

    surfplot(k4x, k5x, rmsez, 'k4', 'k5', 'rmse', 'rmse over k4 and k5 values', False)
    print(rmse([13.68/k1, 186.97/k6], D, S1, XT, q_M, C, C_in,  S1_in, k_hyd))
    return 0;



# define the function for RSE
def rmse(x, D, S1, XT, q_M, C, C_in, S1_in, k_hyd) -> float:
    qClist = []
    for i in range(0, len(D)):
         qClist.append(qCeval(D[i], S1[i], XT[i], q_M[i], C[i], C_in, S1_in, k_hyd, x[0], x[1]))

    return np.sqrt(np.mean(q_C-qClist)**2)

# define the function for the qC
def qCeval(D, S1, XT, q_M, C, C_in, S1_in, k_hyd, CC, DD) -> float:
    return D*(C_in-C) + CC*(D*(S1_in-S1) + k_hyd*XT) + DD*q_M

# utility function for printing and plotting
def fprint(sol, flag, itr, funcall, r2: List):
    print(f"\nSuccess: {flag}\n")
    print("==Results==")
    print(f"Solution: {sol}")
    print(f"k4: {sol[0]*k1}")
    print(f"k5: {sol[1]*k6}")
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
