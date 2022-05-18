# import scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, fsolve
from scipy.optimize import Bounds
from typing import List
from sklearn.metrics import r2_score

from dataimport import*

def minimopt(reg_number,k_hyd, mu_1, mu_2):
    # define the bounds

    # define the initial values ======================== Now one for both
    if reg_number == 4:
        gs = 20
        bd = Bounds([0,0],[1000,1000])
        minsol = minimize(rmse, gs, args = (reg_number, k_hyd, mu_1, mu_2), method = 'SLSQP', tol = 1e-10)
        return minsol.x

    if reg_number == 5:
        gs = 200
        bd = Bounds([0,0],[1000,1000])
        minsol = minimize(rmse, gs, args = (reg_number, k_hyd, mu_1, mu_2), method = 'SLSQP', tol = 1e-10)
        return minsol.x


    if reg_number == 6:
        gs = [1,50]
        bd = Bounds([0,0],[1000,1000])
        minsol = minimize(rmse, gs, args = (reg_number, k_hyd, mu_1, mu_2), method = 'SLSQP', bounds = bd, tol = 1e-10)
        return minsol.x

    if reg_number == 7:
        gs = [1,5]
        bd = Bounds([0,0],[1000,1000])
        minsol = minimize(rmse, gs, args = (reg_number, k_hyd, mu_1, mu_2), method = 'SLSQP', bounds = bd, tol = 1e-10)
        return minsol.x

def rmse(x, reg_number, k_hyd, mu_1, mu_2) -> float:
    Loclist = []

    if reg_number == 4:
        y4exp = D*(S1_in - S1) + k_hyd*XT
        for i in range(0, len(XT)):
                Loclist.append(y4eval(X_1[i], mu_1[i], x[0]))

        return np.sqrt(np.mean(y4exp-Loclist)**2)


    if reg_number == 5:
        y5exp = q_M/X_2
        for i in range(0, len(XT)):
            Loclist.append(y5eval( mu_2[i], x[0]))

        return np.sqrt(np.mean(y5exp-Loclist)**2)


    if reg_number == 6:
        for i in range(0, len(D)):
            Loclist.append(qMeval(D[i], S1[i], S2[i], XT[i], S1_in, S2_in, k_hyd, x[0], x[1]))

        return np.sqrt(np.mean(q_M-Loclist)**2)

    if reg_number == 7:

        for i in range(0, len(D)):
            Loclist.append((qCeval(D[i], S1[i], XT[i], q_M[i], C[i], C_in, S1_in, k_hyd, x[0], x[1])))

        return np.sqrt(np.mean(q_C-Loclist)**2)


def y4eval(X_1, mu_1, k1) -> float:
    return k1*X_1*mu_1

def y5eval(mu_2, k6) -> float:
    return k6*mu_2

def qMeval(D, S1, S2, XT, S1_in, S2_in, k_hyd, AA, BB) -> float:
    return AA*D*(S2_in - S2) + AA*BB*((S1_in-S1) + k_hyd*XT)

def qCeval(D, S1, XT, q_M, C, C_in, S1_in, k_hyd, CC, DD) -> float:
    return D*(C_in-C) + CC*(D*(S1_in-S1) + k_hyd*XT) + DD*q_M
