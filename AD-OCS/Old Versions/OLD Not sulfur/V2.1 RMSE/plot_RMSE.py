
import numpy as np
import matplotlib.pyplot as plt

from Identification_RMSE import*
from dataimport import*

# Plot of k1 regression

print(f'k1: {k1}')
Y_k1 = D*(S1_in - S1) + k_hyd*XT
X_k1 = mu_1*X_1
plt.figure('K1 regression')
plt.plot(X_k1, k1*X_k1, label = 'Fitted')
plt.plot(X_k1, Y_k1, 'o', label = 'Exp. Data')
plt.legend()


# Plot of k6 regression

print(f'k6: {k6}')
Y_k6 = q_M/X_2
X_k6 = mu_2
plt.figure('K6 regression')
plt.plot(X_k6, k6*X_k6, label = 'Fitted')
plt.plot(X_k6, Y_k6, 'o', label = 'Exp. Data')
plt.legend()


# Plot of first ratios regression

plt.figure('AA = k6/k3, BB = k2/k1')
plt.plot(X_k6, k6*X_k6, label = 'Fitted')
plt.plot(X_k6, Y_k6, 'o', label = 'Exp. Data')
plt.legend()
plt.show()