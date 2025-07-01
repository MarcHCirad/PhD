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


def plotTrajectory1D(inputFileNames,
                saveFig = True, saveFileName="",
                xlabel = "", y1label="", y2label = "", y3label="",
                toPlot=False, fontsize = 12, markersize = 20,
                lineLabels = [], titleLegend = ""):

    if saveFileName == "":
        saveFileName = inputFileNames[0][0:-4] + ".png"
    if len(lineLabels) < len(inputFileNames):
        lineLabels = lineLabels + ["" for _ in range(len(inputFileNames)- len(lineLabels))]

    fig, (ax1,ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True)

    for k, inputFileName in enumerate(inputFileNames):

        data = pd.read_csv(inputFileName)
        data = data.to_numpy()
        listParamT = np.float64(data[1:, 0])
        listParamY1 = np.float64(data[1:, 1])
        listParamY2 = np.float64(data[1:, 2])
        listParamY3 = np.float64(data[1:, 3])
        line1, = ax1.plot(listParamT, listParamY1)
        line2, = ax2.plot(listParamT, listParamY2, label=lineLabels[k])
        line3, = ax3.plot(listParamT, listParamY3)


    ax3.set_xlabel(xlabel, fontsize=fontsize)
    ax1.set_ylabel(y1label, fontsize=fontsize)
    ax2.set_ylabel(y2label, fontsize=fontsize)
    ax3.set_ylabel(y3label, fontsize=fontsize)

    ax1.tick_params(labelsize=fontsize)
    ax2.tick_params(labelsize=fontsize)
    ax3.tick_params(labelsize=fontsize)
    ax2.legend(fontsize=fontsize, title = titleLegend, title_fontsize = fontsize)

    if toPlot:
        plt.show()

    if saveFig:
        fig.savefig(saveFileName)

    return fig, (ax1, ax2, ax3)

def plotTrajectory3D(inputFileNames,
                saveFig = True, saveFileName="",
                xlabel = "", ylabel="", zlabel = "",
                toPlot=False, fontsize = 12, markersize = 20,
                lineLabels = [], titleLegend = ""):

    if saveFileName == "":
        saveFileName = inputFileNames[0][0:-4] + ".png"
    if len(lineLabels) < len(inputFileNames):
        lineLabels = lineLabels + ["" for _ in range(len(inputFileNames)- len(lineLabels))]

    fig, ax = plt.subplots(1, subplot_kw={"projection": "3d"})

    for k, inputFileName in enumerate(inputFileNames):

        data = pd.read_csv(inputFileName)
        data = data.to_numpy()
        listParamT = np.float64(data[1:, 0])
        listParamX = np.float64(data[1:, 1])
        listParamY = np.float64(data[1:, 2])
        listParamZ = np.float64(data[1:, 3])
        line, = ax.plot(listParamX, listParamY, listParamZ, label=lineLabels[k])
        ax.scatter(listParamX[0], listParamY[0], listParamZ[0], marker='+', s=markersize, c=line.get_color())
        ax.scatter(listParamX[-1], listParamY[-1], listParamZ[-1], marker='o', s=markersize, c=line.get_color())


    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_zlabel(zlabel, fontsize=fontsize)

    ax.tick_params(labelsize=fontsize)
    ax.legend(fontsize=fontsize, title = titleLegend, title_fontsize = fontsize)

    if toPlot:
        plt.show()

    if saveFig:
        fig.savefig(saveFileName)

    return fig, ax

def plotPhasePlane(inputFileNames,
                saveFig = True, saveFileName="",
                xlabel = "", ylabel="", transformData = lambda x:x,
                toPlot=False, fontsize = 12, markersize = 20,
                lineLabels = [], lineColors=[], titleLegend = ""):
    
    if saveFileName == "":
        saveFileName = inputFileNames[0][0:-4] + ".png"
    if len(lineLabels) < len(inputFileNames):
        lineLabels = lineLabels + ["" for _ in range(len(inputFileNames)- len(lineLabels))]
    if len(lineColors) < len(inputFileNames):
        lineColors = lineColors + ["black" for _ in range(len(inputFileNames)- len(lineColors))]

    fig, ax = plt.subplots(1)
    lines = ["-.","-.",":",":"]
    linecycler = cycle(lines)

    for k, inputFileName in enumerate(inputFileNames):
        data = pd.read_csv(inputFileName)
        data = data.to_numpy()
        data2D, data3D = transformData(data)
        listParamX2D = np.float64(data2D[1:, 0])
        listParamY2D = np.float64(data2D[1:, 1])

        listParamX3D = np.float64(data3D[1:, 0])
        listParamY3D = np.float64(data3D[1:, 1])

        line, = ax.plot(listParamX2D, listParamY2D, label="2D", color=lineColors[k], linestyle=next(linecycler), zorder=0)
        ax.scatter(listParamX2D[0], listParamY2D[0], marker='+', s=markersize, c=line.get_color())
        ax.scatter(listParamX2D[-1], listParamY2D[-1], marker='o', s=4*markersize, c=line.get_color())

        line, = ax.plot(listParamX3D, listParamY3D, label="3D", color=lineColors[k+1], linestyle=next(linecycler), zorder=0)
        ax.scatter(listParamX3D[0], listParamY3D[0], marker='+', s=markersize, c=line.get_color())
        ax.scatter(listParamX3D[-1], listParamY3D[-1], marker='o', s=4*markersize, c=line.get_color())


    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)

    ax.tick_params(labelsize=fontsize)
    ax.legend(fontsize=fontsize, title = titleLegend, title_fontsize = fontsize)

    if toPlot:
        plt.show()

    if saveFig:
        fig.savefig(saveFileName)

    return fig, ax

def separateSimu(tab):
    nrow, _ = tab.shape
    rslt2D = np.zeros((nrow, 2), dtype=float)
    rslt3D = np.zeros((nrow, 2), dtype=float)
    rslt2D = tab[:, 1:3]
    rslt3D = tab[:, 4:6]
            
    return rslt2D, rslt3D

pathWD = str(pathlib.Path(__file__).parent.resolve())
input = pathWD + "/result.csv"

myFig, myAx = plotPhasePlane([input], toPlot=False, saveFig=False,
                                transformData=separateSimu,
                                xlabel = "$H_D$", ylabel="$F_W$",
                               fontsize=20, markersize = 40, lineColors=['black', 'grey'],
                               lineLabels=["$0$"],
                               titleLegend="$\\mathcal{I} = $"
                            )

plt.show(block = True)




