#inertial oscillations exercise #1 ESM
#euler example script

import numpy as np
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
eulerErrArr = np.zeros(numTimeSteps)

#initial conditions
##sim init cond
uArr[0] = 10.0

#matching that for the analytical answer comparison
uAnaArr[0] = uArr[0]


#mathing
##euler forward scheme
###loop over time
for nt in range(1,numTimeSteps):
	timeArr[nt] = timeArr[nt - 1] + dt
	uArr[nt] = uArr[nt-1] + (fCor * vArr[nt-1] * dt)
	vArr[nt] = vArr[nt-1] - (fCor * uArr[nt-1] * dt)

	#getting the analytical answer for this step
	uAnaArr[nt] = uAnaArr[0] * np.cos(fCor * timeArr[nt])

	#comparing the two to see the error
	eulerErrArr[nt] = np.abs(uArr[nt] - uAnaArr[nt])



#plotting
##cleanup crew
plt.close("all")

##plotting the estimated oscillation vs the actual
plt.figure("Oscillation: Euler fwd vs analytical")
plt.title("Oscillation: Euler fwd vs analytical")
plt.xlabel("time\n(s)")
plt.ylabel("u\n(m/s)")
plt.plot(timeArr, uArr, '.b', label="Euler")
plt.plot(timeArr, uAnaArr, "r", label="Analytical")
plt.legend()

#plotting the error (difference between the model and analytical solution)
plt.figure("Error comparison")
plt.title("error")
plt.xlabel("time\n(s)")
plt.ylabel("error (m/s)")
plt.plot(timeArr, eulerErrArr, ".g")

##Don't forget to show 'em if you got 'em
plt.show()

##one for the comparison of model to the analytical solution
##one to show the error growth over time
##perhaps one for energy