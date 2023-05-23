# Author: Alex Rigido

# import libraries & set plot styling
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("ggplot")

# define parameters
x_min = 0       # (m)
x_max = 1       # (m)
time_start = 0  # start time at 0
time_end = 5    # end time
dt = 0.01       # time step (s)
dx = 0.02       # step size (m)
g = 9.81        # gravity

# Time stamps for plotting
time_stamps = np.arange(time_start, time_end + 1)

# set up functions
def wave_height(x, A, mu, sig):
    ''' This function returns the displacement of the wave assuming
        the undisturbed surface of the wave is 0.'''
    n = A * np.exp(-((x - mu) ** 2) / sig ** 2)
    return n

def bottom_depth(x):
    ''' This function return the depth of the bottom surface
        assuming the undisturbed surface of the wave is 0.'''
    if (x =='flat'):
        h = 1/100
    else:
        h = 0.001 + 0.1 * np.e**(-7*x)
    return h


def conservative_SWE(H, U):
    ''' This function is the flux-conservative form of the 1D
        swallow water equations.'''
    # LHS
    u = U[:,0]  # solving for u
    n = U[:,1]  # time index
    # RHS
    F = np.zeros(U.shape)
    F[:,0] = (u**2)/2 + n*g
    F[:,1] = (H+n)*u
    return F

# Define the time & space arrays
x = np.arange(x_min, x_max + dx, dx)
t = np.arange(time_start, time_end + dt, dt)
A = 0.002
mu = 0.5
sig = 0.05
# main algorithm (2-Step Lax-Wendroff)
u = np.reshape(x*0, (x.shape[0], 1))
n = np.reshape(wave_height(x, A, mu, sig), (x.shape[0], 1))
U = np.concatenate((u, n), axis=1)
x_values = []

for n in range(0, len(t)):
    # complete full step
    F = conservative_SWE(bottom_depth('flat'), U)
    U_half_step = np.zeros((U.shape[0] - 1, U.shape[1]))
    U[0] = np.array([0, U[0, 1] - (dt/dx)*(F[1, 1] - F[0, 1])])
    U[-1] = np.array([0, U[-1, 1] - (dt/dx)*(F[-1, 1] - F[-2, 1])])
    # complete half step
    U_half_step[0] = 0.5 * (U[1] + U[0]) - (dt / dx) * (F[1] - F[0])
    U_half_step[-1] = 0.5 * (U[-2] + U[-1]) - (dt / dx) * (F[-1] - F[-2])
    U_half_step[1:-1] = 0.5 * ((U[2:-1] + U[1:-2]) - (dt / dx) * (F[2:-1] - F[1:-2]))
    F_half_step = conservative_SWE(bottom_depth('flat'), U_half_step)
    # Store values
    U[1:-1] = U[1:-1] - (dt/dx)*(F_half_step[1:] - F_half_step[:-1])
    x_values.append(np.copy(U[:, 1]))

# Plot results at time stamps.
plt.figure()
plt.title("1-D Shallow Waves (Flat Bottom Topography)")
plt.xlabel("Wave Position (m)")
plt.ylabel("Wave Height (m)")
for stamp in time_stamps:
    plt.plot(x, x_values[np.where(t == stamp)[0][0]], lw="2", label=(str(stamp) + " Seconds"))
plt.legend()
# plt.savefig("flat_bottom_SWE.png")
plt.show()
