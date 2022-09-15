import numpy as np 

X2 = np.linspace(1,2,20)

tspan = np.linspace(0, 10, 20)

t = 1
A = X2[tspan<t]
# get the corresponding index of t in tspan
idx = np.where(tspan>t)
print(idx)  # get the last index

