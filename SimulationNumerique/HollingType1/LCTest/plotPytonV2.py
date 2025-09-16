import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pathlib

pathWD = str(pathlib.Path(__file__).parent.resolve())


inputFile = pathWD + "/resultLC" + ".csv"

beginLC = 0
dt = 0.01
endLC = int(20000 / dt)
beginLC = int(0/dt)
endLC = int(20000/dt)

plotLCExtInt = False

if plotLCExtInt:
        inputFileExt = pathWD + "/resultLCExt" + ".csv"
        inputFileInt = pathWD + "/resultLCInt" + ".csv"
        inputs = [inputFile, inputFileExt, inputFileInt]
else:
     inputs = [inputFile]


plot = ['orbit', "time"]
# plot = ['orbit']  ## pour avoir uniquement les orbites
# plot = ['time'] ## pour avoir uniquement le LC en fonction du temps

LC = 0
Ext = 1
Int = 2
colorLC = 'red'
colorInt = 'orange'
colorExt = 'blue'

fontsize = 30
markersize = 300
linewidth = 1
linewidthLC = 3


datasT = []
datasFW = []
datasHD = []
datasHW = []

for file in inputs:
    data = pd.read_csv(file)
    data = data.to_numpy()
    data = data
    nrow, _ = data.shape

    dataT = np.zeros((nrow, 1), dtype=float)
    data3D = np.zeros((nrow, 2), dtype=float)
    dataT = data[:, 0]
    data3D = data[:, 4:7]
    HDVal3D = data3D[:, 0]
    FWVal3D = data3D[:, 1]
    HWVal3D = data3D[:, 2]

    datasT.append(dataT)
    datasFW.append(FWVal3D)
    datasHD.append(HDVal3D)
    datasHW.append(HWVal3D)

print("min FW LC : ", np.min(datasFW[LC]))
print("min FW LC : ", np.max(datasFW[LC]))
last = LC
print(datasHD[last][-1], datasFW[last][-1], datasHW[last][-1])


if "orbit" in plot:
    fig, ax = plt.subplots(1)
    line, = ax.plot(datasHD[LC][beginLC:endLC]+ datasHW[LC][beginLC:endLC], datasFW[LC][beginLC:endLC], color= colorLC, linestyle='--', zorder=2, linewidth= linewidthLC,
                    label = 'LC')
    ax.scatter(datasHD[LC][0] + datasHW[LC][0], datasFW[LC][0], marker='d', s = markersize, c = colorLC,
            )
    if plotLCExtInt:
        line, = ax.plot(datasHD[Int] + datasHW[Int], datasFW[Int], color= colorInt, linestyle='--', zorder=0, linewidth=linewidth,
                        label = 'Starts inside the LC')
        line, = ax.plot(datasHD[Ext] + datasHW[Ext], datasFW[Ext], color = colorExt, linestyle = '--', zorder=1, linewidth=linewidth,
                        label = 'Starts outside the LC')
        ax.scatter(datasHD[Int][0] + datasHW[Int][0], datasFW[Int][0], marker='d', s = markersize, c = colorInt,
            )
        ax.scatter(datasHD[Ext][0] + datasHW[Ext][0], datasFW[Ext][0], marker='d', s = markersize, c = colorExt,
            )
        
    ax.set_xlabel('$H_D + H_W$', fontsize=fontsize)
    ax.set_ylabel('$F_W$', fontsize=fontsize)

    ax.tick_params(labelsize=fontsize)
    legend = ax.legend(fontsize=fontsize, title = '', title_fontsize = fontsize, loc = 'upper right')
    for line in legend.get_lines():
        line.set_linewidth(linewidthLC)
    plt.show()

if "time" in plot:
    fig1, ax1 = plt.subplots(nrows=1, ncols=1, sharex=True)

    ax1.plot(datasT[LC][0:endLC-beginLC], datasHD[LC][beginLC:endLC] + datasHW[LC][beginLC:endLC], label="$H_D + H_W$", color='black',
            linewidth = linewidthLC,
            linestyle = ':')
    ax1.plot(datasT[LC][0:endLC-beginLC], datasFW[LC][beginLC:endLC], label="$F_W$", color='black',
            linewidth = linewidthLC,
            linestyle = '--')

    ax1.set_xlabel('$t$', fontsize=fontsize)
    ax1.tick_params(labelsize=fontsize)
    ax1.legend(fontsize=fontsize)
    ax1.set_xlim(0 * dt, (endLC-beginLC) * dt)


    plt.show()