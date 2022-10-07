import numpy as np
import pprint as pp
from matplotlib import pyplot as plt
# Level dynamics function
def level_dynamics(t, D, Qin, SR, h0, t_change):
  t_loc = t-t_change 
  return h0*np.exp(-D*t_loc) + (1-np.exp(-D*t_loc))*Qin/(D*SR)

# Available data
t_span = np.linspace(0,100,100) # [d] Time span
Qin = 170           # m3/d  - volumetric inlet flowrate
t_change = [20, 60]   # [d]   - time of change in the inlet flowrate - N changes
scale = [1, 1.8, 1]   # [-]   - scale factor for the inlet flowrate  - N+1 factors (the first is always one)
index = 0

D = 0.05              # 1/d   - Dilution rate
h0 = 10.823               # m     - initial level in the reactor
DR = 20               # m     - reactor diameter
SR = np.pi/4 * DR**2  # m2    - reactor surface
hmax = 15
hmin = 5

# Simulation solution
t_change_loc = t_span[0]
h = np.empty(len(t_span))
for i in range(len(t_span)):
  t = t_span[i]   
  if t < t_change[index]:     
    Qin_loc = Qin*scale[index]  

  elif t >= t_change[-1]:
    Qin_loc = Qin*scale[-1]
    t_change_loc = t_change[-1] 
    h0 = h[len(t_span[t_span < t_change_loc])-1]

  else:       
    t_change_loc = t_change[index] 
    index = index + 1
    Qin_loc = Qin*scale[index]
    h0 = h[i-1]
 
  h[i] = level_dynamics(t, D, Qin_loc, SR, h0, t_change_loc)
  if h[i] > hmax:
      print('*** Kittemmuort è troppo pieno ***')
      input('Press Enter to continue')
  elif h[i] < hmin:
      print('*** ittemmuort è troppo vuoto ***')
      input('Press Enter to continue')
  print('time: ', "{:.0f}".format(t_span[i]), ' d - level: ', "{:.3f}".format(h[i]), ' m - flow: ', "{:.1f}".format(Qin_loc), ' m3/d' , 'h0: ', "{:.3f}".format(h0))

# Visualization
plt.figure()
plt.plot(t_span, h, 'b', linewidth = 2)
#plt.plot(t_span, [hmax for _ in range(len(h))], 'k--')
plt.xlabel('time [d]')
plt.ylabel('level [m]')
plt.grid(True)
plt.show()