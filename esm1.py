#inertial oscillations exercise #1 ESM
#euler example script

import numpy as np
import matplotlib.pyplot as plt

#consts
numTimeSteps = 800
dt = 200.0 #s?
fCor = 0.001 #Hz

#data struct setup
##time step array and analytical solution array
timeArr = np.zeros(numTimeSteps)
uAnaArr = np.zeros(numTimeSteps)

##arrays for euler fwd method
uEulerArr = np.zeros(numTimeSteps)
vEulerArr = np.zeros(numTimeSteps)
eulerErrArr = np.zeros(numTimeSteps)

##arrays for the leapfrog method
uLeapFrogArr = np.zeros(numTimeSteps)
vLeapFrogArr = np.zeros(numTimeSteps)
leapFrogErrArr = np.zeros(numTimeSteps)

##arrays for the Heun method
uHeunArr = np.zeros(numTimeSteps)
vHeunArr = np.zeros(numTimeSteps)
heunErrArr = np.zeros(numTimeSteps)

##arrays for the Matsuno method
uMatsunoArr = np.zeros(numTimeSteps)
vMatsunoArr = np.zeros(numTimeSteps)
matsunoErrArr = np.zeros(numTimeSteps)



#initial conditions
##sim init cond are all the same
uEulerArr[0] = 10.0
uLeapFrogArr[0] = uEulerArr[0]
uHeunArr[0] = uEulerArr[0]
uMatsunoArr[0] = uEulerArr[0]

#matching analytical answer comparison
uAnaArr[0] = uEulerArr[0]


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
	#
	# Heun scheme
	uHeunPred = uHeunArr[nt-1] + (fCor * vHeunArr[nt-1] * dt)
	vHeunPred = vHeunArr[nt-1] - (fCor * uHeunArr[nt-1] * dt)
	uHeunArr[nt] = uHeunArr[nt-1] + (fCor * (vHeunPred + vHeunArr[nt-1]) * 0.5 * dt)
	vHeunArr[nt] = vHeunArr[nt-1] - (fCor * (uHeunPred + uHeunArr[nt-1]) * 0.5 * dt)
	heunErrArr[nt] = np.abs(uHeunArr[nt] - uAnaArr[nt])

	#matsuno
	# Matsuno scheme
	if nt == 1:
	    vMatsunoArr[1] = vEulerArr[nt]
	    print("Matsuno IC #2 trig")
	else:
	    uMatsunoArr[nt] = uMatsunoArr[nt-2] + 2 * (fCor * vMatsunoArr[nt-1] * dt)
	    vMatsunoArr[nt] = vMatsunoArr[nt-1] - (fCor * (uMatsunoArr[nt-2] + uMatsunoArr[nt]) * dt)
	matsunoErrArr[nt] = np.abs(uMatsunoArr[nt] - uAnaArr[nt])


	#getting the analytical answer for this step
	uAnaArr[nt] = uAnaArr[0] * np.cos(fCor * timeArr[nt])

	#comparing the schemes to see the error
	eulerErrArr[nt] = np.abs(uEulerArr[nt] - uAnaArr[nt])
	leapFrogErrArr[nt] = np.abs(uLeapFrogArr[nt] - uAnaArr[nt])


#plotting
##cleanup crew
plt.close("all")

##plotting the estimated oscillation vs the actual
plt.figure("Oscillation: Euler fwd vs analytical")
plt.title("Euler fwd vs analytical")
plt.xlabel("time\n(s)")
plt.ylabel("u\n(m/s)")
plt.plot(timeArr, uEulerArr, '.b', label="Euler")
plt.plot(timeArr, uAnaArr, "r", label="Analytical")
# !!! put the other methods here too
plt.legend()

##plotting the other schemes
plt.figure("Oscillation: others vs analytical")
plt.title("Cd vs Heun vs Matsuno vs analytical")
plt.xlabel("time\n(s)")
plt.ylabel("u\n(m/s)")
plt.plot(timeArr, uLeapFrogArr, ".c", label="Leapfrog")
plt.plot(timeArr, uHeunArr, ".m", label="Heun")
#plt.plot(timeArr, uMatsunoArr, ".y", label="Matsuno")
plt.plot(timeArr, uAnaArr, "k", label="Analytical")
plt.legend()

##plot heun scheme


##plot matsuno scheme

#plotting the error (difference between the model and analytical solution)
plt.figure("Error comparison")
plt.title("Euler error")
plt.xlabel("time\n(s)")
plt.ylabel("error (m/s)")
plt.plot(timeArr, eulerErrArr, ".r", label="Euler")
plt.legend()

plt.figure("other errors")
plt.title("Cd vs Heun vs Matsuno error")
plt.xlabel("time\n(s)")
plt.ylabel("error (m/s)")
plt.plot(timeArr, leapFrogErrArr, ".g", label="leapfrog")
plt.plot(timeArr, heunErrArr, ".r", label="Heun")
#plt.plot(timeArr, matsunoErrArr, ".b", label="Matsuno")
plt.legend()


##Don't forget to show 'em if you got 'em
plt.show()

##one to show the error growth over time
##perhaps one for energy