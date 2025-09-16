import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pathlib

pathWD = str(pathlib.Path(__file__).parent.resolve())

inputFileLC = pathWD + "/resultLC" + ".csv"
inputFileLCBeta = pathWD + "/resultLCBeta" + ".csv"
inputFileLCTheta005 = pathWD + "/resultLCTheta005" + ".csv"
inputFileLCTheta01 = pathWD + "/resultLCTheta01" + ".csv"
inputFileLCTheta01B = pathWD + "/resultLCTheta01B" + ".csv"
inputFileLCTheta1 = pathWD + "/resultLCTheta1" + ".csv"

plot = ["orbit"]#, "time"]

LC = 0
LCBeta = 1
LC1 = 1
LC2 = 2
LC2B = 3
LC3 = 4

inputs3D = [inputFileLC]#, inputFileLCBeta]
inputs2D = [inputFileLCTheta005, inputFileLCTheta01, inputFileLCTheta01B, inputFileLCTheta1]

colorLC = 'red'
colorLC1 = 'orange'
colorLC2 = 'blue'
colorLC2B = 'black'
colorLC3 = 'gray'

fontsize = 30
markersize = 300
linewidth = 1
linewidthLC = 3
begin = 0
dt = 0.01
end = int(20000 / dt)
beginLC = int(0/dt)
endLC = int(200/dt)

datasT = []
datasFW = []
datasHD = []
datasHW = []

for file in inputs3D:
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

for file in inputs2D:
    data = pd.read_csv(file)
    data = data.to_numpy()
    data = data
    nrow, _ = data.shape

    dataT = np.zeros((nrow, 1), dtype=float)
    data2D = np.zeros((nrow, 2), dtype=float)
    dataT = data[:, 0]
    data2D = data[:, 4:7]
    HDVal2D = data2D[:, 0]
    FWVal2D = data2D[:, 1]
    HWVal2D = data2D[:, 2]

    datasT.append(dataT)
    datasFW.append(FWVal2D)
    datasHD.append(HDVal2D)
    datasHW.append(HWVal2D)

print("min FW LC : ", np.min(datasFW[LC]))
print("min FW LC : ", np.max(datasFW[LC]))
last = LC2
print(datasHD[last][-1], datasFW[last][-1], datasHW[last][-1])


if "orbit" in plot:
    fig, ax = plt.subplots(1)
    line, = ax.plot(datasHD[LC][beginLC:endLC], datasFW[LC][beginLC:endLC], color= colorLC, linestyle='--', zorder=2, linewidth= linewidthLC,
                    label = '$\\theta=0, \\lambda_F = 0.06$')
    # line, = ax.plot(datasHD[LCBeta][beginLC:endLC], datasFW[LCBeta][beginLC:endLC], linestyle='--', zorder=2, linewidth= linewidthLC,
    #                 label = '$\\theta=0, \\lambda_F = 0.06, \\beta > 0$')
    line, = ax.plot(datasHD[LC1][beginLC:endLC], datasFW[LC1][beginLC:endLC], color= colorLC1, linestyle='--', zorder=2, linewidth= linewidthLC,
                    label = '$\\theta=0.005, \\lambda_F = 0.06$')
    line, = ax.plot(datasHD[LC2][beginLC:endLC], datasFW[LC2][beginLC:endLC], color= colorLC2, linestyle='--', zorder=2, linewidth= linewidthLC,
                    label = '$\\theta=0.1, \\lambda_F = 0.06$')
    line, = ax.plot(datasHD[LC2B][beginLC:endLC], datasFW[LC2B][beginLC:endLC], color= colorLC2B, linestyle='--', zorder=2, linewidth= linewidthLC,
                    label = '$\\theta=0.1, \\lambda_F = 0.003$')
    line, = ax.plot(datasHD[LC3][beginLC:endLC], datasFW[LC3][beginLC:endLC], color= colorLC3, linestyle='--', zorder=2, linewidth= linewidthLC,
                    label = '$\\theta=1, \\lambda_F = 0.0006$')
    # ax.scatter(datasHD[Int][0], datasFW[Int][0], marker='d', s = markersize, c = colorInt,
    #         )
    # ax.scatter(datasHD[Ext][0] + datasHW[Ext][0], datasFW[Ext][0], marker='d', s = markersize, c = colorExt,
    #         )

    ax.set_xlabel('$H_D$', fontsize=fontsize)
    ax.set_ylabel('$F_W$', fontsize=fontsize)

    ax.tick_params(labelsize=fontsize)
    legend = ax.legend(fontsize=fontsize, title = '', title_fontsize = fontsize, loc = 'upper right')
    for line in legend.get_lines():
        line.set_linewidth(linewidthLC)
    plt.show()

if "time" in plot:
    fig1, ax1 = plt.subplots(nrows=1, ncols=1, sharex=True)
    dt = 0.01
    begin = int(19000 / dt)
    end = int(20000/dt)

    plotting = LCBeta

    ax1.plot(datasT[plotting][0:end-begin], datasHD[plotting][begin:end], label="$H_D + H_W$", color='black',
            linewidth = linewidthLC,
            linestyle = ':')
    ax1.plot(datasT[plotting][0:end-begin], datasFW[plotting][begin:end], label="$F_W$", color='black',
            linewidth = linewidthLC,
            linestyle = '--')

    ax1.set_xlabel('$t$', fontsize=fontsize)
    ax1.tick_params(labelsize=fontsize)
    ax1.legend(fontsize=fontsize)
    ax1.set_xlim(0 * dt, (end-begin) * dt)


    plt.show()