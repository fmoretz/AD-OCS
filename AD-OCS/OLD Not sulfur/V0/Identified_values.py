# IDENTIFICATION
import numpy as np

# Kinetics
C_d = [0.1, 0.1]
mu_max = [0.315740612955049, 0.132093656852006]        # [1/d]      Max biomass growth rate
Ks     = [0.337542594773213, 3.790224063348594]        # [g/L]      Half saturation constant
KI2    = abs(193.756138)              # [mmol/L]   Inhibition constant for S2
kd  = np.multiply(C_d,mu_max)

# Physical
kLa = 31.15                      # [1/d]      L/G transfer rate

# Yield coefficients: [CODdeg, VFAprof, VFAcons, CO2prod(1), CO2prod(2), CH4prod, hydrolysis]
k = [18.5102608361823,	802.438017101711,	958.808000569365,	5.01810987833723,	200.647370296197,	252.899306266204,	5.00241137750437]                # [ --  mmol/g mmol/g ....] ==================================================
