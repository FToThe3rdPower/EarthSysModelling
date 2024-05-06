#inertial oscillations exercise #1 ESM
#euler example script

#!!! put in a func so I can run it over different time lengths
#!!! and vary the dt

import numpy as np
import matplotlib.pyplot as plt

#sim vars & consts
numTimeSteps = 800
dt = 200.0 #s?
fCor = 0.001 #Hz
initCond  = 10.0

#data struct setup
##time step array and analytical solution array
timeArr = np.zeros(numTimeSteps)
uAnaArr = np.zeros(numTimeSteps)

##arrays for euler fwd method
uEulerArr = np.zeros(numTimeSteps)
vEulerArr = np.zeros(numTimeSteps)
eulerErrArr = np.zeros(numTimeSteps)
eulerEngArr = np.zeros(numTimeSteps)

##arrays for the leapfrog method
uLeapFrogArr = np.zeros(numTimeSteps)
vLeapFrogArr = np.zeros(numTimeSteps)
leapFrogErrArr = np.zeros(numTimeSteps)
leapFrogEngArr = np.zeros(numTimeSteps)

##arrays for the Heun method
uHeunArr = np.zeros(numTimeSteps)
vHeunArr = np.zeros(numTimeSteps)
heunErrArr = np.zeros(numTimeSteps)
heunEngArr = np.zeros(numTimeSteps)


##arrays for the Matsuno method
uMatsunoArr = np.zeros(numTimeSteps)
vMatsunoArr = np.zeros(numTimeSteps)
matsunoErrArr = np.zeros(numTimeSteps)
matsunoEngArr = np.zeros(numTimeSteps)


#initial conditions
##sim init cond are all the same
uEulerArr[0] = initCond
uLeapFrogArr[0] = initCond
uHeunArr[0] = initCond
uMatsunoArr[0] = initCond

#matching analytical answer comparison
uAnaArr[0] = initCond


#mathing
##loop over time
for nt in range(1,numTimeSteps):
	#euler forward scheme
	timeArr[nt] = timeArr[nt - 1] + dt
	uEulerArr[nt] = uEulerArr[nt-1] + (fCor * vEulerArr[nt-1] * dt)
	vEulerArr[nt] = vEulerArr[nt-1] - (fCor * uEulerArr[nt-1] * dt)

	#leapfrog scheme
	##give it the second point from Euler
	if nt == 1:
		uLeapFrogArr[1] = uEulerArr[nt]
		vLeapFrogArr[1] = vEulerArr[nt]

	##proceed with the rest of the lf scheme
	else:
		uLeapFrogArr[nt] = fCor*vLeapFrogArr[nt-1]*2*dt + uLeapFrogArr[nt-2] 
		vLeapFrogArr[nt] = -fCor*uLeapFrogArr[nt-1]*2*dt + vLeapFrogArr[nt-2]
	
	#heun scheme
	uHeunPred = uHeunArr[nt-1] + (fCor * dt * vHeunArr[nt-1])
	vHeunPred = vHeunArr[nt-1] - (fCor * dt * uHeunArr[nt-1])
	uHeunArr[nt] = uHeunArr[nt-1] + (fCor * 0.5 * dt * (vHeunPred + vHeunArr[nt-1]))
	vHeunArr[nt] = vHeunArr[nt-1] - (fCor * 0.5 * dt * (uHeunPred + uHeunArr[nt-1]))

	#matsuno
	uMatsunoPred = uMatsunoArr[nt-1] + (fCor * dt * vMatsunoArr[nt-1])
	vMatsunoPred = vMatsunoArr[nt-1] - (fCor * dt * uMatsunoArr[nt-1])
	uMatsunoArr[nt] = uMatsunoArr[nt-1] + (fCor * dt * vMatsunoPred)
	vMatsunoArr[nt] = vMatsunoArr[nt-1] - (fCor * dt * uMatsunoPred)


	#calculating energies
	eulerEngArr[nt] = 0.5*np.abs(uEulerArr[nt]**2 + vEulerArr[nt]**2)
	leapFrogEngArr[nt] = 0.5*np.abs(uLeapFrogArr[nt]**2 + vLeapFrogArr[nt]**2)
	heunEngArr[nt] = 0.5*np.abs(uHeunArr[nt]**2 + vHeunArr[nt]**2)
	matsunoEngArr[nt] = 0.5*np.abs(uMatsunoArr[nt]**2 + vMatsunoArr[nt]**2)

	#getting the analytical answer for this step
	uAnaArr[nt] = uAnaArr[0] * np.cos(fCor * timeArr[nt])

	#comparing the schemes to see the error
	eulerErrArr[nt] = np.abs(uEulerArr[nt] - uAnaArr[nt])
	leapFrogErrArr[nt] = np.abs(uLeapFrogArr[nt] - uAnaArr[nt])
	heunErrArr[nt] = np.abs(uHeunArr[nt] - uAnaArr[nt])
	matsunoErrArr[nt] = np.abs(uMatsunoArr[nt] - uAnaArr[nt])


#plotting
##cleanup crew
plt.close("all")

##plotting the estimated oscillation vs the actual
plt.figure("Oscillation: Euler fwd vs analytical")
plt.title("Euler fwd vs Matsuno with analytical")
plt.xlabel("time\n(s)")
plt.ylabel("u\n(m/s)")
plt.plot(timeArr, uEulerArr, '.b', label="Euler")
plt.plot(timeArr, uAnaArr, "r", label="Analytical")
plt.legend()

##plotting the other schemes
plt.figure("Oscillation: others vs analytical")
plt.title("leapfrog, Heun, matsuno vs analytical")
plt.xlabel("time\n(s)")
plt.ylabel("u\n(m/s)")
plt.plot(timeArr, uLeapFrogArr, ".c", label="Leapfrog")
plt.plot(timeArr, uHeunArr, ".m", label="Heun")
plt.plot(timeArr, uMatsunoArr, ".y", label="Matsuno")
plt.plot(timeArr, uAnaArr, "k", label="Analytical")
plt.legend()

##energy plots
plt.figure("Oscillations: energy")
plt.title("Euler fwd vs cd vs Heun vs Matsuno energy conservation")
plt.xlabel("time\n(s)")
plt.ylabel("Energy")
plt.semilogy(timeArr, eulerEngArr, label="Euler")
plt.semilogy(timeArr, leapFrogEngArr, label="Leapfrog")
plt.semilogy(timeArr, heunEngArr, label="Heun")
plt.semilogy(timeArr, matsunoEngArr, label="Matsuno")
plt.legend()

##plotting the error (difference between the model and analytical solution)
plt.figure("Error comparison")
plt.title("Euler  error")
plt.xlabel("time\n(s)")
plt.ylabel("error (m/s)")
plt.semilogy(timeArr, eulerErrArr, ".r", label="Euler")
plt.legend()

plt.figure("Leapfrog Heun Matsuno vs analytical errors")
plt.title("Leapfrog, Heun, Matsuno, vs analytical error")
plt.xlabel("time\n(s)")
plt.ylabel("error (m/s)")
plt.plot(timeArr, leapFrogErrArr, ".g", label="leapfrog")
plt.plot(timeArr, heunErrArr, ".r", label="Heun")
plt.plot(timeArr, matsunoErrArr, ".b", label="Matsuno")
plt.legend()


##Don't forget to show 'em if you got 'em
plt.show()

##one to show the error growth over time
##perhaps one for energy