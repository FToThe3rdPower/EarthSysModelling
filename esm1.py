#inertial oscillations exercise #1 ESM
#euler example script

import numpy as np
import matplotlib.pyplot as plt

#consts
numTimeSteps = 720
dt = 300.0 #s?
fCor = 0.00011 #Hz

#data struct setup
##time step array and analytical solution array
timeArr = np.zeros(numTimeSteps)
uAnaArr = np.zeros(numTimeSteps)

##arrays for euler fwd method
uEulerArr = np.zeros(numTimeSteps)
vEulerArr = np.zeros(numTimeSteps)
eulerErrArr = np.zeros(numTimeSteps)

##arrays for the Heun method
uHeunArr = np.zeros(numTimeSteps)
vHeunArr = np.zeros(numTimeSteps)
heunErrArr = np.zeros(numTimeSteps)

##arrays for the Matsuno method
uMatsunoArr = np.zeros(numTimeSteps)
vMatsunoArr = np.zeros(numTimeSteps)
matsunoErrArr = np.zeros(numTimeSteps)



#initial conditions
##sim init cond ar all the same
uEulerArr[0] = 10.0
uHeunArr[0] = uEulerArr[0]
uMatsunoArr[0] = uEulerArr[0]

#matching analytical answer comparison
uAnaArr[0] = uEulerArr[0]


#mathing
##loop over time
for nt in range(1,numTimeSteps):
	#euler forward scheme
	timeArr[nt] = timeArr[nt - 1] + dt
	uEulerArr[nt] = uEulerArr[nt-1] + (fCor * vEulerArr[nt-1] * dt
	vEulerArr[nt] = vEulerArr[nt-1] - (fCor * uEulerArr[nt-1] * dt)

	#getting the analytical answer for this step
	uAnaArr[nt] = uAnaArr[0] * np.cos(fCor * timeArr[nt])

	#comparing the two to see the error
	eulerErrArr[nt] = np.abs(uEulerArr[nt] - uAnaArr[nt])


	#heun scheme
	#!!!



#plotting
##cleanup crew
plt.close("all")

##plotting the estimated oscillation vs the actual
plt.figure("Oscillation: Euler fwd vs analytical")
plt.title("Oscillation: Euler fwd vs analytical")
plt.xlabel("time\n(s)")
plt.ylabel("u\n(m/s)")
plt.plot(timeArr, uEulerArr, '.b', label="Euler")
plt.plot(timeArr, uAnaArr, "r", label="Analytical")
plt.legend()

##plot heun scheme

##plot matsuno scheme

#plotting the error (difference between the model and analytical solution)
plt.figure("Error comparison")
plt.title("error")
plt.xlabel("time\n(s)")
plt.ylabel("error (m/s)")
plt.plot(timeArr, eulerErrArr, ".g", label="Euler")

##Don't forget to show 'em if you got 'em
plt.show()

##one for the comparison of model to the analytical solution
##one to show the error growth over time
##perhaps one for energy