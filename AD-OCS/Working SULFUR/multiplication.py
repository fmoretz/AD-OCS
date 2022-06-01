''' Perform an array multiplication and show the result'''
import numpy as np 
from Influent import*
from dataimport import*
from deviationscopy2 import deviation_check as dc
scale = [1, 1, 1, 1, 1]
scale2 = [1, 1, 1, 1.2, 1]

dev = y_in_0*scale


dev2 = y_in_0*scale2


res = dc(35,y_in_0)
print(res)