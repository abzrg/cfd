import warnings
warnings.filterwarnings("ignore")

import sys

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def main():
    data = read_data()

    nx = len(data)
    x = np.linspace(0.0, 2.0 * np.pi, nx)

    df = pd.DataFrame({"x": x, "y": data})

    plt.figure(figsize=(10, 6))
    plt.plot(df["x"], df["y"], marker="o", markersize=3, linestyle="-", color="k")
    plt.title(r"Initial Condition: $u_0$")
    plt.xlabel(r"$x$")
    plt.ylabel(r"$u$")
    plt.grid(True)
    plt.savefig("U0.png")
    plt.show()


def usage():
    print("Usage: plot1d <field>")


def read_data():
    fpath: str
    if len(sys.argv) == 2:
        fpath = sys.argv[1]
    else:
        usage()
        raise SystemExit(-1)

    data = pd.read_csv(fpath, header=None).squeeze().tolist()
    return data


main()
