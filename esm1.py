#inertial oscillations exercise #1 ESM
#euler example script

import numpy as numpy
import matplotlib.pyplot as plt

#consts
numTimeSteps = 720
dt = 300.0 #s?
fCor = 0.00011 #Hz

#data struct setup
timeArr = np.zeros(numTimeSteps)
uAnaArr = np.zeros(numTimeSteps)
uArr = np.zeros(numTimeSteps)
vArr = np.zeros(numTimeSteps)


#initial conditions
##sim init cond
u[0] = 10.0

#matching that for the analytical answer comparison
uAna[0] = u[0]



#mathing: euler forward scheme


#plotting
##one for the comparison of model to the analytical solution
##one to show the error growth over time
##perhaps one for energy