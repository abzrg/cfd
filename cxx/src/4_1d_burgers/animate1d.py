import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
import re


def extract_number(filename):
    """Extract time from 'time_step_TIME.txt' file names."""
    match = re.search(r"timestep_(\d+\.?\d+)\.txt", filename)
    return float(match.group(1)) if match else 0


def read_data(file):
    """Read data from a timestep result file ('time_step_TIME.txt')."""
    data = pd.read_csv(file, header=None).squeeze().tolist()
    return data


# List all files and sort them numerically
file_list = [
    f for f in os.listdir(".") if f.startswith("timestep_") and f.endswith(".txt")
]
file_list.sort(key=extract_number)  # Sort files based on the numerical part

# Create lists to store data and times
data_list = []
times = []

for filename in file_list:
    time = extract_number(filename)
    times.append(time)
    data_list.append(read_data(filename))

# Create x-domain
n = len(data_list[0]) if data_list else 0
x = np.linspace(0.0, 2.0*np.pi, n)

# Create the figure and axis for plotting
fig, ax = plt.subplots(figsize=(10, 6))
(line,) = ax.plot(x, data_list[0], marker="o", linestyle="-")


# Show the animation
ax.set_xlabel("x")
ax.set_ylabel("u")
plt.grid(True)


# Function to update the plot
def update(frame):
    line.set_ydata(data_list[frame])
    # Not sure why this does not get updated in the backend window,
    # but gets updated in the saved animation (gif or mp4)
    ax.set_title(f"Time: {times[frame]}")
    return (line,)


# Create the animation
ani = animation.FuncAnimation(
    fig, update, frames=len(data_list), blit=True, repeat=True
)
ani.save("animation.gif")
# ani.save("animation.mp4")

plt.show()
