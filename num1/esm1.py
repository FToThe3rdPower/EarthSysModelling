#inertial oscillations exercise #1 ESM
#euler example script

#!!! put in a func so I can run it over different time lengths
#!!! and vary the dt

import numpy as np
import matplotlib.pyplot as plt

#sim vars & consts
numberOfTimeSteps = 500
dt = 200 #s
frequency = 0.00011 #Hz
initCond  = 10.0
sizeOfDTStep = 100
dtStart = sizeOfDTStep
dtEnd = 1000

rangeOfDt = range(dtStart, dtEnd, sizeOfDTStep)
sizeOfDTsweep = (dtEnd - dtStart)/sizeOfDTStep

#data struct setup
##time step array and analytical solution array
timeArr = np.zeros(numberOfTimeSteps)
uAnaArr = np.zeros(numberOfTimeSteps)
vAnaArr = np.zeros(numberOfTimeSteps)
engAnaArr = np.zeros(numberOfTimeSteps)

##arrays for euler fwd method
uEulerArr = np.zeros(numberOfTimeSteps)
vEulerArr = np.zeros(numberOfTimeSteps)
eulerErrArr = np.zeros(numberOfTimeSteps)
eulerEngArr = np.zeros(numberOfTimeSteps)

##arrays for the leapfrog method
uLeapFrogArr = np.zeros(numberOfTimeSteps)
vLeapFrogArr = np.zeros(numberOfTimeSteps)
leapFrogErrArr = np.zeros(numberOfTimeSteps)
leapFrogEngArr = np.zeros(numberOfTimeSteps)

##arrays for the Heun method
uHeunArr = np.zeros(numberOfTimeSteps)
vHeunArr = np.zeros(numberOfTimeSteps)
heunErrArr = np.zeros(numberOfTimeSteps)
heunEngArr = np.zeros(numberOfTimeSteps)

##arrays for the Matsuno method
uMatsunoArr = np.zeros(numberOfTimeSteps)
vMatsunoArr = np.zeros(numberOfTimeSteps)
matsunoErrArr = np.zeros(numberOfTimeSteps)
matsunoEngArr = np.zeros(numberOfTimeSteps)


#arrays for sepecific method err, will be avg'd
eulerErrVdt = np.empty(sizeOfDTsweep)
leapFrogErrVdt = np.empty(sizeOfDTsweep)
heunErrVdt = np.empty(sizeOfDTsweep)
matsunoErrVdt = np.empty(sizeOfDTsweep)


#initial conditions
##sim init cond are all the same
uEulerArr[0] = initCond
uLeapFrogArr[0] = initCond
uHeunArr[0] = initCond
uMatsunoArr[0] = initCond

#matching analytical answer comparison
uAnaArr[0] = initCond


#mathing
##this function calculates the oscillation analytically
##along with estimation based on different schemes and their errors
#to spit out an array of arrays with the speficif info to be plotted
def FUNCytown(numTimeSteps, dt, fCor):
	##loop over time
	for nt in range(1,numTimeSteps):
		timeArr[nt] = timeArr[nt - 1] + dt
		
		#getting the analytical answer for this step
		uAnaArr[nt] = uAnaArr[0] * np.cos(fCor * timeArr[nt])
		vAnaArr[nt] = vAnaArr[0] * np.sin(fCor * timeArr[nt])
		
		#euler forward scheme
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
		engAnaArr[nt] = 0.5*(initCond**2)
		eulerEngArr[nt] = 0.5*np.abs(uEulerArr[nt]**2 + vEulerArr[nt]**2)
		leapFrogEngArr[nt] = 0.5*np.abs(uLeapFrogArr[nt]**2 + vLeapFrogArr[nt]**2)
		heunEngArr[nt] = 0.5*np.abs(uHeunArr[nt]**2 + vHeunArr[nt]**2)
		matsunoEngArr[nt] = 0.5*np.abs(uMatsunoArr[nt]**2 + vMatsunoArr[nt]**2)


		#comparing the schemes to see the error
		eulerErrArr[nt] = np.abs(uEulerArr[nt] - uAnaArr[nt])
		leapFrogErrArr[nt] = np.abs(uLeapFrogArr[nt] - uAnaArr[nt])
		heunErrArr[nt] = np.abs(uHeunArr[nt] - uAnaArr[nt])
		matsunoErrArr[nt] = np.abs(uMatsunoArr[nt] - uAnaArr[nt])

	#return an array of the arrays we'll need so we can run it multiple times to get info about varying the dt
	return np.array([timeArr, uAnaArr, uEulerArr, uLeapFrogArr, uHeunArr, uMatsunoArr, eulerEngArr, leapFrogEngArr, heunEngArr, matsunoEngArr, eulerErrArr, leapFrogErrArr, heunErrArr, matsunoErrArr])

#running
arr4funcAtDT = FUNCytown(numberOfTimeSteps, dt, frequency)

#plotting
##cleanup crew
plt.close("all")

#!!! add function for plotting


##plotting the estimated oscillation vs the actual
plt.figure("Zonal Velocity v time: Leapfrog, Heun, Matsuno vs analytical")
plt.title("Zonal Velocity v time:\nLeapfrog, Heun, Matsuno vs analytical\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("time\n(s)")
plt.ylabel("u\n(m/s)")
plt.plot(arr4funcAtDT[0], arr4funcAtDT[1], "k", label="Analytical")
plt.plot(arr4funcAtDT[0], arr4funcAtDT[2], '.b', label="Euler")
plt.plot(arr4funcAtDT[0], arr4funcAtDT[3], "xc", label="Leapfrog")
plt.plot(arr4funcAtDT[0], arr4funcAtDT[4], ".m", label="Heun")
plt.plot(arr4funcAtDT[0], arr4funcAtDT[5], ".y", label="Matsuno")
plt.legend()

##energy plots
plt.figure("Kinetic energy v time: Leapfrog v Heun")
plt.title("Kinetic Energy v time:\nLeapfrog vs Heun\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("time\n(s)")
plt.ylabel("Energy")
plt.plot(arr4funcAtDT[0], leapFrogEngArr, "xr", label="Leapfrog")
plt.plot(arr4funcAtDT[0], heunEngArr, ".b", label="Heun")
plt.plot(arr4funcAtDT[0], engAnaArr, "k", label="Analytic")

plt.legend()

plt.figure("Kinetic Energy v time: Euler")
plt.title("Kinetic Energy v time: Euler\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("time\n(m/s)")
plt.ylabel("Energy")
plt.semilogy(arr4funcAtDT[0], eulerEngArr,".g", label="Euler")
plt.semilogy(arr4funcAtDT[0], engAnaArr, "k", label="Analytic")
plt.legend()

plt.figure("Kinetic Energy v time: Matsuno")
plt.title("Kinetic Energy v time: Matsuno\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("time\n(m/s)")
plt.ylabel("Energy")
plt.plot(arr4funcAtDT[0], matsunoEngArr,".b", label="Matsuno")
plt.plot(arr4funcAtDT[0], engAnaArr, "k", label="Analytic")

plt.legend()

##plotting the error (difference between the model and analytical solution)
plt.figure("Zonal Velocity error v time: Euler method error")
plt.title("Euler Velocity error\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("time\n(s)")
plt.ylabel("error (m/s)")
plt.semilogy(arr4funcAtDT[0], eulerErrArr, ".r", label="Euler")
plt.legend()

plt.figure("Zonal Velocity error v time: Leapfrog, Heun vs analytical")
plt.title("Zonal Velocity error v time:\nLeapfrog, Heun vs analytical\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("time\n(s)")
plt.ylabel("error (m/s)")
plt.plot(arr4funcAtDT[0], leapFrogErrArr, "xg", label="leapfrog")
plt.plot(arr4funcAtDT[0], heunErrArr, ".r", label="Heun")
plt.legend()

plt.figure("Zonal Velocity error v time: Matsuno vs analytical")
plt.title("Zonal Velocity error v time:\nMatsuno vs analytical\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("time\n(s)")
plt.ylabel("error (m/s)")
plt.plot(arr4funcAtDT[0], matsunoErrArr, ".b", label="Matsuno")
plt.legend()

#loop over a few dt options and 
for nDT in rangeOfDt:
	#get the avg error as a func of nDT for each method 
	arrToChop = FUNCytown(numberOfTimeSteps, nDT, frequency)
	
	#grabbing a method's avg error array out of the arrArr
	eulerErrVdt[(nDT/sizeOfDTStep)-1] = np.mean(arrToChop[10])
	leapFrogErrVdt[(nDT/sizeOfDTStep-1)] = np.mean(arrToChop[11])
	heunErrVdt[(nDT/sizeOfDTStep)-1] = np.mean(arrToChop[12])
	matsunoErrVdt[(nDT/sizeOfDTStep)-1] = np.mean(arrToChop[13])

plt.figure("avg error v dt: euler")
plt.title("Model error v Timestep\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("dt\n(x100 s)")
plt.ylabel("np.avg(Model's error array)")
plt.semilogy(eulerErrVdt, ".", label="euler")
plt.legend()

plt.figure("avg error v dt: heun")
plt.title("Model error v Timestep\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("dt\n(x100 s)")
plt.ylabel("np.avg(Model's error array)")
plt.semilogy(heunErrVdt, "xg", label="heun")
plt.legend()

plt.figure("avg error v dt: leapfrog and matsuno")
plt.title("Model error v Timestep\n#Tsteps = "
	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
plt.xlabel("dt\n(x100 s)")
plt.ylabel("np.avg(Model's error array)")
plt.plot(leapFrogErrVdt, "xr", label="leapfrog")
plt.plot(matsunoErrVdt, ".b", label="matsuno")
plt.legend()

#plt.figure("Err v dt: mat")
#plt.title("Model error v Timestep\n#Tsteps = "
#	+ str(numberOfTimeSteps) + " dt = " + str(dt) + " fCor = " + str(frequency) + " IC= " + str(initCond))
#plt.xlabel("dt\n(x100)")
#plt.ylabel("np.avg(Model's error array)")
#plt.plot(matsunoErrVdt, "x", label="matsuno")
#plt.legend()

##Don't forget to show 'em if you got 'em
plt.show()