# Author: Alex Rigido
# Please note, this will be the vectorized simulation for a
# faster implementation and better run time.

# import libraries & set plot styling
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("ggplot")

# define parameters
x_min = 0       # (m)
x_max = 1      # (m)
time_start = 0  # start time at 0
time_end = 5    # end time
dt = 0.001      # time step (s)
dx = 1/150      # step size (m)
g = 9.81        # gravity

# Time stamps for plotting
time_stamps = np.arange(time_start, time_end + 1)

# Define the position, speed, & space arrays
x = np.arange(x_min, x_max + dx,dx)
u = np.zeros_like(x)
n = np.zeros_like(x)

# define the bottom_topography for simulation
height = np.clip(0.0005 + (0.01*(1-np.tanh(8*np.pi*(x-0.5)))), 0, 2)
height[:] = 0.1*(np.exp(-x*7)) + 0.001

# set up functions
def wave_height(x):
    ''' This function returns the displacement of the wave assuming
        the undisturbed surface of the wave is 0.'''
    A = 0.0005  # m
    mu = 0.0    # m
    sig = 0.1   # m
    n = A * np.e **(-((x-mu)**2)/sig**2)
    return n

def flux_step(U, gravity, height):
    ''' This function is the flux step of the 2-Step Lax-Wendroff method and
        returns flux calculated at a given time-step, where input U = U(u,n).'''
    flux = np.zeros_like(U)
    flux[0] = gravity * U[1] + 0.5*U[0]**2
    flux[1] = (U[1] + height)*U[0]
    return flux

def lax_step(vector, dx, dt, gravity, height):
    ''' This function is the lax step of the 2-Step Lax-Wendroff method and
            returns the calculated lax-step for an inputted (2 X N) vector.'''
    # compute midpoint averages
    vector_midpoint = 0.5*(vector[:, :-1] + vector[:, 1:])
    height_midpoint = 0.5*(height[1:] + height[:-1])

    # compute flux on boundaries
    boundary_flux = flux_step(vector, gravity, height)

    # compute centered differences
    flux_midpoint = (boundary_flux[:,1:] - boundary_flux[:,:-1])
    vector_half = vector_midpoint - dt * flux_midpoint / (2*dx)

    # compute flux at half time & space
    flux_half_height = flux_step(vector_half, gravity, height_midpoint)
    flux_midpoint_center = (flux_half_height[:,1:] - flux_half_height[:,:-1])
    flux_midpoint_step = np.zeros_like(vector)
    flux_midpoint_step[:,1:-1] = flux_midpoint_center

    # set boundary conditions
    flux_midpoint_step[0,0] = 0
    flux_midpoint_step[0,-1] = 0
    flux_midpoint_step[1,0] = (boundary_flux[1,1] - boundary_flux[1,0])
    flux_midpoint_step[1,-1] = (boundary_flux[1,-1] - boundary_flux[1,-2])

    return vector - dt * flux_midpoint_step / dx

# main algorithm (2-Step Lax-Wendroff)
n = wave_height(x)      # initial wave height
u[:] = 0                # assumed initial wave speed 0
U = np.vstack([u, n])   # initial state of the system U(u,n)

# define the time array for iteration and initialize storage
time_array = np.arange(time_start, time_end, dt)
result = np.zeros((time_array.size, U.shape[0], U.shape[1]))

# iteration through time array
for i,t in enumerate(time_array):
    U = lax_step(U, dx, dt, g, height)
    result[i] = U

# visualize simulation
plt.figure()
plt.title("Tsunami Simulation")
plt.xlabel("Wave Position (m)")
# add the wave height data
plt.ylabel("Wave Height (m)")
plt.ylim(-0.001,0.001)
for stamp in time_stamps:
    plt.plot(x, result[int(stamp // dt), 1], lw=2, label=(str(stamp) + " Seconds"))
plt.legend(loc=4)
# add the bottom topography
plt.twinx()
plt.plot(x,-height, color='black', label = "Bottom Profile")
plt.ylabel("Bottom Height (m)")
plt.ylim(-0.1,0.1)
plt.legend()
# plt.savefig("tsunami_simulation_SWE.png")
plt.show()