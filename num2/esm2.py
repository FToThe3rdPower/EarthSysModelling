#second numerics exercise
import numpy as np
import matplotlib.pyplot as plt

#vars
##Sim choices
schemeFlag = "leapFrog"

u0 = 0.01 #10/1000 # km/s but 5km/min seems pretty fast, that's 300km/hr...
c0 = 1
L = 2500 #km

#Default dx = 25km, dckdt overflow occurs for dt > 1000s
dx = 25 #km
#for the loop over dx's to see the effect on the schemes
rangeMax = 500
dxRange = np.arange(100, rangeMax, 100, dtype=int)

dt = 100 #s default 300
dtRangeMax = 500
dtRange = np.arange(100, rangeMax, 100, dtype=int)

ridx = 1375
cidx = 1250
lidx = 1125


##Number of dt, 1 cycle
nt_len =  int((L/u0)/dt)
sigma = u0 * dt/dx #$\leq$ 1

##array setup
t_hour = np.linspace(0, dt/3600 * (nt_len), nt_len + 1)
t_axis = np.linspace(0, dt * nt_len, nt_len + 1)
x_axis = np.linspace(0 + dx/2, L-dx/2, int(L/dx))

u0_tracer_c = np.mod(cidx + u0 * t_axis, L)
u0_tracer_l = np.mod(lidx + u0 * t_axis, L)
u0_tracer_r = np.mod(ridx + u0 * t_axis, L)


#FUNCy town
##Diffusion and dispersion
def EulerUpwind(dx,dt,nt_len):
    nx_len = int(L/dx)
    c_sol = np.zeros((nt_len + 1,nx_len))
    pulse_st_idx = int(np.floor(lidx/dx + 1)-1)
    pulse_ed_idx = int(np.ceil(ridx/dx)-1)
    c_sol[0,pulse_st_idx:pulse_ed_idx + 1] = c0 * np.ones((1,pulse_ed_idx-pulse_st_idx + 1))
    for n in range(1,nt_len + 1):
        c_sol[n,:] = c_sol[n-1,:] - sigma * (c_sol[n-1,:]-np.roll(c_sol[n-1,:],1))
    return c_sol

##diffusion
def LaxWendroff(dx,dt,nt_len):
    nx_len = int(L/dx)
    c_sol = np.zeros((nt_len + 1,nx_len))
    pulse_st_idx = int(np.floor(lidx/dx + 1)-1)
    pulse_ed_idx = int(np.ceil(ridx/dx)-1)
    c_sol[0,pulse_st_idx:pulse_ed_idx + 1] = c0 * np.ones((1,pulse_ed_idx-pulse_st_idx + 1))
    for n in range(1,nt_len + 1):
        c_sol[n,:] = c_sol[n-1,:] - sigma/2 * (np.roll(c_sol[n-1,:],-1)-np.roll(c_sol[n-1,:],1)) + np.power(sigma,2)/2 * (np.roll(c_sol[n-1,:],-1)-2 * c_sol[n-1,:] + np.roll(c_sol[n-1,:],1))
    return c_sol

#spectral derivative
def dckdt(cks0):
    ks = np.arange(0, np.size(cks0))
    return -u0 * cks0 * 1j * 2 * np.pi * ks / L

def Forward(vec0, func, dt, type):
	#local vars in your area lookin to compute
    dvecdt0 = func(vec0)
    vec1_ = vec0 + dvecdt0 * dt
    vec1 =  vec1_

    #specifying scheme flag
    if type != 'leapFrog':
        dvecdt1_ = func(vec1_)
        if type == 'Matsuno':
            vec1 = vec0 + (0 * dvecdt0 + 1 * dvecdt1_) * dt
        elif type == 'Heun':
            vec1 = vec0 + (0.5 * dvecdt0 + 0.5 * dvecdt1_) * dt
    return vec1


#spectral with numpy's fft alg
def Spectral(dx, dt, nt_len):
	nx_len = int(L / dx)
	nj_len = nx_len # // 2 + 1  # Include Nyquist component

	# Arrays for the transform
	c_sol = np.zeros((nt_len + 1, nx_len))
	ck_sol = np.zeros((nt_len + 1, nj_len), dtype=complex)

	#c_sol init stuff
	pulse_st_idx = int(np.floor(lidx/dx + 1)-1)
	pulse_ed_idx = int(np.ceil(ridx/dx)-1)
	c_sol[0, pulse_st_idx:pulse_ed_idx + 1] = c0 * np.ones((1,pulse_ed_idx-pulse_st_idx + 1))


	# Wavenumber
	k = np.fft.fftfreq(nx_len) * (2 * np.pi / L)

	# Initial condition in spectral space
	ck_sol[0, :] = np.fft.fft(c_sol[0, :]) * dx / (nj_len - 1)

	for n in range(1, nt_len + 1):
		# Time evolution in spectral space using Forward func
		ck_sol[n, :] = Forward(ck_sol[n-1, :], dckdt, dt, schemeFlag)

		# Back to real space
		c_sol[n, :] = np.real(np.fft.ifft(ck_sol[n, :] * (nj_len - 1) / dx))

	return c_sol




#calls
c_sol_euler = EulerUpwind(dx,dt,nt_len)
c_sol_LaxWendroff = LaxWendroff(dx,dt,nt_len)



# #plots
# ##euler concentration plat
# plt.figure("eulerC")
# plt.title("euler concentration, dx = "+str(dx))
# plt.plot(x_axis, c_sol_euler[0,:], "o-r", label="0")
# plt.plot(x_axis, c_sol_euler[100,:], "o-g", label="100")
# plt.plot(x_axis, c_sol_euler[200,:], "o-b", label="200")
# plt.plot(x_axis, c_sol_euler[500,:], "o-c", label="500")
# plt.plot(x_axis, c_sol_euler[1000,:], "o-m", label="1000")
# plt.plot(x_axis, c_sol_euler[1500,:], "o-y", label="1500")
# plt.plot(x_axis, c_sol_euler[2000,:], "o-k", label="2000")
# plt.xlabel('x\n(km)')
# plt.ylabel('$C$\n($C_{0}$)')
# plt.legend(loc="upper left")

# ##laxWen con plot
# plt.figure("laxWenC")
# plt.title("laxWen concen. dx = "+str(dx))
# plt.plot(x_axis, c_sol_LaxWendroff[0,:], "o-r", label="0")
# plt.plot(x_axis, c_sol_LaxWendroff[10,:], "o-g", label="10")
# plt.plot(x_axis, c_sol_LaxWendroff[50,:], "o-b", label="50")
# plt.plot(x_axis, c_sol_LaxWendroff[200,:], "o-c", label="200")
# plt.plot(x_axis, c_sol_LaxWendroff[500,:], "o-m", label="500")
# plt.plot(x_axis, c_sol_LaxWendroff[1500,:], "o-y", label="1500")
# plt.plot(x_axis, c_sol_LaxWendroff[2000,:], "o-k", label="2000")
# plt.xlabel('x\n(km)')
# plt.ylabel('$C$\n($C_{0}$)')
# plt.legend()


# #Euler Hovmoller diagram
# plt.figure("eulerHov")
# plt.title("Hovmoller diagram, dx = " + str(dx) +"\nEuler method forward in time and space")

# #normalize the colormap
# cmin = np.amin(c_sol_euler)
# cmax = np.amax(c_sol_euler)

# #make the grid for the trace lines
# xx, tt = np.meshgrid(x_axis, t_axis)

# #plot the results
# colorMesh  = plt.pcolormesh(xx, tt, c_sol_euler, cmap="ocean", vmin=cmin, vmax=cmax)

# #add the trace lines
# tracerline1_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c > cidx)[0]], t_axis[np.where(u0_tracer_c > cidx)[0]], color='r', lw=3)
# tracerline2_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c < cidx)[0]], t_axis[np.where(u0_tracer_c < cidx)[0]], color='r', lw=3)

# colorBar = plt.colorbar(colorMesh)
# colorBar.set_label(label='\n$C$\n($C_{0}$)')
# plt.xlabel('x\n(km)')
# plt.ylabel('time\n(s)')



# #Lax-Wendroff Hovmoller diagram
# plt.figure("laxWendHov")
# plt.title('Hovmoller diagram, dx = ' +str(dx)+' Lax-Wendroff')

# #normalize the colormap
# cmin = np.amin(c_sol_LaxWendroff)
# cmax = np.amax(c_sol_LaxWendroff)

# xx, tt = np.meshgrid(x_axis, t_axis)
# laxColorMesh = plt.pcolormesh(xx, tt, c_sol_LaxWendroff, cmap="ocean",  vmin=cmin, vmax=cmax)
# tracerline1_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c > cidx)[0]], t_axis[np.where(u0_tracer_c > cidx)[0]], color='r', lw=3)
# tracerline2_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c < cidx)[0]], t_axis[np.where(u0_tracer_c < cidx)[0]], color='r', lw=3)

# #colorBar setup
# laxColorBarf = plt.colorbar(laxColorMesh)
# laxColorBarf.set_label('\n$C$\n($C_{0}$)')
# plt.xlabel('x\n(km)')
# plt.ylabel('time\n(s)')



# # Plotting separate figures for each scheme
# plt.figure("eulerDX")

# ## Loop over dx values
# for dx in dxRange:

# 	print(dx)
# 	# Calculate sigma based on the current dx
# 	sigma = u0 * dt / dx

# 	# Calculate solution for Euler 
# 	c_sol_euler = EulerUpwind(dx, dt, nt_len)

# 	# Plot Euler solutions
# 	plt.plot(c_sol_euler[-1, :], label=dx)
# plt.xlabel("dx (km)")
# plt.ylabel("Final Concentration")
# plt.title("Concentration, Euler FWD " + schemeFlag)
# plt.legend()



# # Plotting separate figures for each scheme
# plt.figure("LaxWendroffDX")

# ## Loop over dx values
# for dx in dxRange:

# 	print(dx)
# 	# Calculate sigma based on the current dx
# 	sigma = u0 * dt / dx

# 	# Calculate solution for Euler 
# 	c_sol_laxw = LaxWendroff(dx, dt, nt_len)


# 	# Plot final concentration
# 	plt.plot(c_sol_laxw[-1, :], label=dx)
# plt.xlabel("dx (km)")
# plt.ylabel("Final Concentration")
# plt.title("Concentration, LaxWendroff")
# plt.legend()




# Plotting separate figures for each scheme
plt.figure("eulerDT")

## Loop over dx values
for dt in dtRange:
	nt_len =  int((L/u0)/dt)

	# Calculate sigma based on the current dx
	sigma = u0 * dt / dx

	# Calculate solution for Euler 
	c_sol_euler = EulerUpwind(dx, dt, nt_len)

	# Plot Euler solutions
	plt.plot(c_sol_euler[-1, :], label=dt)
plt.xlabel("dt (s)")
plt.ylabel("Final Concentration")
plt.title("Concentration, Euler FWD " + schemeFlag)
plt.legend()



# Plotting separate figures for each scheme
plt.figure("LaxWendroffDT")

## Loop over dx values
for dt in dtRange:
	nt_len =  int((L/u0)/dt)

	# Calculate sigma based on the current dx
	sigma = u0 * dt / dx

	# Calculate solution for Euler 
	c_sol_laxw = LaxWendroff(dx, dt, nt_len)


	# Plot final concentration
	plt.plot(c_sol_laxw[-1, :], label=dt)
plt.xlabel("dt (s)")
plt.ylabel("Final Concentration")
plt.title("Concentration, LaxWendroff")
plt.legend()


#spectral stuff brushed under the rug
#print(Spectral(dx,dt,nt_len))
#c_sol_Spectral = Spectral(dx,dt,nt_len)

# ##Spectral con plot
# plt.figure("SpectralC " + schemeFlag)
# plt.title("SpectralC " + schemeFlag)
# plt.plot(x_axis, c_sol_Spectral[0,:], "o-r", label="0")
# plt.plot(x_axis, c_sol_Spectral[10,:], "o-b", label="10")
# plt.plot(x_axis, c_sol_Spectral[20,:], "o-g", label="20")
# plt.plot(x_axis, c_sol_Spectral[50,:], "o-c", label="50")
# plt.plot(x_axis, c_sol_Spectral[100,:], "o-m", label="100")
# plt.xlabel('x\n(km)')
# plt.ylabel('$C$\n($C_{0}$)')
# plt.legend()

# #Spectral Hovmoller diagram
# plt.figure("SpectralHov_" + schemeFlag)
# plt.title('Hovmoller diagram, Spectral, ' + schemeFlag)

# #normalize the colormap
# cmin = np.amin(c_sol_Spectral)
# cmax = np.amax(c_sol_Spectral)

# xx, tt = np.meshgrid(x_axis, t_axis)
# colourMesh = plt.pcolormesh(xx, tt, c_sol_Spectral, cmap = "ocean", vmin=cmin, vmax=cmax)
# tracerline1_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c > cidx)[0]], t_axis[np.where(u0_tracer_c > cidx)[0]], color='w', lw=3)
# tracerline2_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c < cidx)[0]], t_axis[np.where(u0_tracer_c < cidx)[0]], color='w', lw=3)

# #color setup
# cb = plt.colorbar(colourMesh)
# cb.set_label(label='\n$C$\n($C_{0}$)')
# plt.xlabel('x\n(km)')
# plt.ylabel('time\n(s)')


#render plots
plt.show()