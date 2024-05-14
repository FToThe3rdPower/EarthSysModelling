#second numerics exercise
import numpy as np
import matplotlib.pyplot as plt

#vars
##Sim choices
schemeFlag = "leapFrog"

u0 = 0.01 #10/1000 # km/s but 5km/min seems pretty fast, that's 300km/hr...
c0 = 1
L = 2500 #km

#D#efault dx = 25km, dckdt overflow occurs for dt > 1000s
dx = 25 #km
dt = 100 #s 300

##Number of dt, 1 cycle
nt_len =  int((L/u0)/dt)

##array setup
t_hour = np.linspace(0, dt/3600 * (nt_len), nt_len + 1)
t_axis = np.linspace(0, dt * nt_len, nt_len + 1)
x_axis = np.linspace(0 + dx/2, L-dx/2, int(L/dx))
cidx = 1250
u0_tracer_c = np.mod(cidx + u0 * t_axis, L)
lidx = 1125
u0_tracer_l = np.mod(lidx + u0 * t_axis, L)
ridx = 1375
u0_tracer_r = np.mod(ridx + u0 * t_axis, L)


#FUNCy town
##Diffusion and dispersion
def EulerUpwind(dx,dt,nt_len):
    sigma = u0 * dt/dx #$\leq$ 1
    nx_len = int(2500/dx)
    c_sol = np.zeros((nt_len + 1,nx_len))
    pulse_st_idx = int(np.floor(1125/dx + 1)-1)
    pulse_ed_idx = int(np.ceil(1375/dx)-1)
    c_sol[0,pulse_st_idx:pulse_ed_idx + 1] = c0 * np.ones((1,pulse_ed_idx-pulse_st_idx + 1))
    for n in range(1,nt_len + 1):
        c_sol[n,:] = c_sol[n-1,:] - sigma * (c_sol[n-1,:]-np.roll(c_sol[n-1,:],1))
    return c_sol

##diffusion
def LaxWendroff(dx,dt,nt_len):
    sigma = u0 * dt/dx #$\leq$ 1
    nx_len = int(2500/dx)
    c_sol = np.zeros((nt_len + 1,nx_len))
    pulse_st_idx = int(np.floor(1125/dx + 1)-1)
    pulse_ed_idx = int(np.ceil(1375/dx)-1)
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

def Spectral(dx,dt,nt_len):
    sigma = u0 * dt/dx
    nx_len = int(2500/dx)
    nj_len = int(nx_len/2) + 1
    c_sol = np.zeros((nt_len + 1,nx_len))

    #spectral coeff
    ck_sol = 1j *  np.zeros((nt_len + 1,nj_len))
    pulse_st_idx = int(np.floor(1125/dx + 1)-1)
    pulse_ed_idx = int(np.ceil(1375/dx)-1)
    c_sol[0,pulse_st_idx:pulse_ed_idx + 1] = c0 * np.ones((1,pulse_ed_idx-pulse_st_idx + 1))
    cks0 = 1j * np.ones((1,nj_len))
    cks0[0,0] = np.mean(c_sol[0,:])
    
    jvec = np.reshape(np.arange(1, nx_len + 1),(nx_len,1))
    kvec = np.reshape(np.arange(1, nj_len),(1,nj_len-1))
    ft_coef = (2/nx_len) * np.exp(-1j * (2 * np.pi/nx_len) * np.matrix(jvec * kvec))
    cks0[0,1:] = np.matrix(c_sol[0,:] * ft_coef)
    ck_sol[0,:] = cks0
    
    for n in range(1,nt_len + 1):
        ck_sol[n,:] = Forward(ck_sol[n-1,:], dckdt, dt, schemeFlag)

    jvec1 = np.reshape(np.arange(1, nx_len + 1),(1,nx_len))
    kvec1 = np.reshape(np.arange(0, nj_len),(nj_len,1))
    c_sol = np.real(np.matrix(ck_sol * (np.exp( 1j * (2 * np.pi/nx_len) * np.matrix(kvec1 * jvec1)))))
    return c_sol



#calls
c_sol_euler = EulerUpwind(dx,dt,nt_len)
c_sol_LaxWendroff = LaxWendroff(dx,dt,nt_len)
c_sol_Spectral = Spectral(dx,dt,nt_len)
print(c_sol_Spectral)



#plots
plt.figure("euler")

#euler plot
plt.title("euler 0")
plt.plot(x_axis, c_sol_euler[0,:], ".r")
plt.xlabel('x\n(km)')
plt.ylabel('$C$\n($C_{0}$)')



#Euler Hovmoller diagram
plt.figure("eulerHov")
plt.title('Hovmoller diagram\nEuler method forward in time and space')

#normalize the colormap
cmin = np.amin(c_sol_euler)
cmax = np.amax(c_sol_euler)

#make the grid for the trace lines
xx, tt = np.meshgrid(x_axis, t_axis)

#plot the results
colorMesh  = plt.pcolormesh(xx, tt, c_sol_euler, cmap="ocean", vmin=cmin, vmax=cmax)

#add the trace lines
tracerline1_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c > cidx)[0]], t_axis[np.where(u0_tracer_c > cidx)[0]], color='r', lw=3)
tracerline2_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c < cidx)[0]], t_axis[np.where(u0_tracer_c < cidx)[0]], color='r', lw=3)

colorBar = plt.colorbar(colorMesh, orientation="vertical")
colorBar.set_label(label='\n$C$\n($C_{0}$)')
plt.xlabel('x\n(km)')
plt.ylabel('time\n(s)')



#Lax-Wendroff Hovmoller diagram
plt.figure("laxWendHov")
plt.title('Hovmoller diagram, Lax-Wendroff')

#normalize the colormap
cmin = np.amin(c_sol_LaxWendroff)
cmax = np.amax(c_sol_LaxWendroff)

xx, tt = np.meshgrid(x_axis, t_axis)
laxColorMesh = plt.pcolormesh(xx, tt, c_sol_LaxWendroff, cmap="ocean",  vmin=cmin, vmax=cmax)
tracerline1_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c > cidx)[0]], t_axis[np.where(u0_tracer_c > cidx)[0]], color='r', lw=3)
tracerline2_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c < cidx)[0]], t_axis[np.where(u0_tracer_c < cidx)[0]], color='r', lw=3)

#colorBar setup
laxColorBarf = plt.colorbar(laxColorMesh, orientation="vertical")
laxColorBarf.set_label('\n$C$\n($C_{0}$)')
plt.xlabel('x\n(km)')
plt.ylabel('time\n(s)')



#Spectral Hovmoller diagram
plt.figure("SpectralHov")
plt.title('Hovmoller diagram, Spectral, ')

#normalize the colormap
cmin = np.amin(c_sol_Spectral)
cmax = np.amax(c_sol_Spectral)

xx, tt = np.meshgrid(x_axis, t_axis)
colourMesh = plt.pcolormesh(xx, tt, c_sol_Spectral, cmap = "ocean", vmin=cmin, vmax=cmax)
tracerline1_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c > cidx)[0]], t_axis[np.where(u0_tracer_c > cidx)[0]], color='w', lw=3)
tracerline2_c = plt.plot(u0_tracer_c[np.where(u0_tracer_c < cidx)[0]], t_axis[np.where(u0_tracer_c < cidx)[0]], color='w', lw=3)

#color setup
cb = plt.colorbar(colourMesh, orientation="vertical")
cb.set_label(label='\n$C$\n($C_{0}$)')
plt.xlabel('x\n(km)')
plt.ylabel('time\n(s)')



#render plots
plt.show()