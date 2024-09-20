import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob

def read_data(filename):
    return np.loadtxt(filename)

file_list = sorted(glob.glob("timestep_*.txt"))

initial_data = read_data(file_list[0])
nx, ny = initial_data.shape

fig, ax = plt.subplots()
cax = ax.imshow(initial_data, cmap='coolwarm', origin='lower')
fig.colorbar(cax)

nx = 81
dx = 2.0 / (nx - 1)
sigma = 0.2
dt = sigma * dx
def update(frame):
    data = read_data(file_list[frame])
    cax.set_array(data)
    ax.set_title(f'Time = {frame * dt:.2f}')
    return cax,

ani = FuncAnimation(fig, update, frames=len(file_list), interval=100, blit=True)

plt.show()
