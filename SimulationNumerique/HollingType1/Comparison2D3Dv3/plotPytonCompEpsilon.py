import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib import ticker
from itertools import cycle
import pathlib

pathWD = str(pathlib.Path(__file__).parent.resolve())

epsilon = 10
inputFile = pathWD + "/resultEps" + str(epsilon) + ".csv"

legend = ''
if epsilon == 10:
    legend = "10^{-1}"
    color = 'blue'
    linestyle = '-.'
elif epsilon == 100:
    legend = "10^{-2}"
    color='grey'
    linestyle = '-.'
elif epsilon == 365:
    legend = "\\dfrac{1}{365}"
    color = 'teal'
    linestyle = '-.'
elif epsilon == 1000 or epsilon == 10000:
    legend = "10^{-3}" if epsilon == 1000 else "10^{-4}"
    color ='green'
    linestyle = '-.'

HDEq = 124.52818
FWEq = 23.14

data = pd.read_csv(inputFile)
data = data.to_numpy()
data = data
nrow, _ = data.shape

dataT = np.zeros((nrow, 1), dtype=float)
data2D = np.zeros((nrow, 2), dtype=float)
data3D = np.zeros((nrow, 2), dtype=float)
dataT = data[:, 0]
data2D = data[:, 1:3]
data3D = data[:, 4:6]

HDVal2D = data2D[:, 0]
FWVal2D = data2D[:, 1]

HDVal3D = data3D[:, 0]
FWVal3D = data3D[:, 1]

fig, ax = plt.subplots(1)
fontsize = 30
markersize = 300
linewidth = 5

line, = ax.plot(HDVal2D, FWVal2D, color='orange', linestyle='--', zorder=0, linewidth=linewidth, label="2D system")
# line, = ax.plot(HDVal3D, FWVal3D, label="3D system, $\\varepsilon = " + legend + "$", color = color, linestyle = linestyle, zorder=1, linewidth=linewidth)

ax.scatter(HDVal3D[0], FWVal3D[0], marker='d', s = markersize, c = 'red', label='Initial Condition')
ax.scatter(HDEq, FWEq, marker='o', s = markersize, color='red', label='$EE^{HF_W}$')

ax.set_xlabel('$H_D$', fontsize=fontsize)
ax.set_ylabel('$F_W$', fontsize=fontsize)

ax.tick_params(labelsize=fontsize)
ax.legend(fontsize=fontsize, title = '', title_fontsize = fontsize, loc = 'upper right')

plt.show()
