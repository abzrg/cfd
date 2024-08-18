import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data from the file
data = pd.read_csv("final.txt", header=None).squeeze().tolist()

# Number of data points
nx = len(data)

# Create an x-domain from 0.0 to 2.0
x = np.linspace(0.0, 2.0, nx)

df = pd.DataFrame({"x": x, "y": data})

plt.figure(figsize=(10, 6))
plt.plot(df["x"], df["y"], marker="o", linestyle="-", color="b")
plt.title(r"$u(x,\,t)\quad t = \text{endTime}$")
plt.xlabel(r"$x$")
plt.ylabel(r"$u$")
plt.grid(True)
plt.show()
