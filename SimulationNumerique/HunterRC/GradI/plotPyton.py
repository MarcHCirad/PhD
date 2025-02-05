import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib import ticker
from itertools import cycle


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
    lines = ["-","-.",":",":"]
    linecycler = cycle(lines)

    for k, inputFileName in enumerate(inputFileNames):

        data = pd.read_csv(inputFileName)
        data = data.to_numpy()
        dataTransformed = transformData(data)
        listParamX = np.float64(dataTransformed[1:, 0])
        listParamY = np.float64(dataTransformed[1:, 1])

        line, = ax.plot(listParamX, listParamY, label=lineLabels[k], color=lineColors[k], linestyle=next(linecycler), zorder=0)
        ax.scatter(listParamX[0], listParamY[0], marker='+', s=markersize, c=line.get_color())
        ax.scatter(listParamX[-1], listParamY[-1], marker='o', s=4*markersize, c=line.get_color())


    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)

    ax.tick_params(labelsize=fontsize)
    ax.legend(fontsize=fontsize, title = titleLegend, title_fontsize = fontsize)

    if toPlot:
        plt.show()

    if saveFig:
        fig.savefig(saveFileName)

    return fig, ax

  
def plotContinuousBifurcationFile(inputFileName,
                saveFig = True, saveFileName="",
                xlabel = "", ylabel="", colorbarTitle = "",
                normName = "",
                toPlot=False, fontsize = 12):

    if saveFileName == "":
        saveFileName = inputFileName[0:-4] + ".png"

    data = pd.read_csv(inputFileName)
    data = data.to_numpy()
    listParamX = np.float64(data[0, 1:])
    listParamY = np.float64(data[1:, 0])
    values = np.array(data[1:, 1:], dtype=float)

    fig, ax = plt.subplots(1)

    mNorm = None
    match normName:
        case "SymLogNorm":
            mNorm = colors.SymLogNorm(linthresh = 0.0001, vmin=values.min(), vmax=values.max())
        case "PowerNorm":
            mNorm = colors.PowerNorm(gamma=0.05)
        case "LogNorm":
            mNorm=colors.LogNorm(vmin=values.min(), vmax=values.max())
        case _:
            mNorm = None

    im = ax.pcolormesh(listParamX, listParamY, values, norm = mNorm)
    cbar = fig.colorbar(im)
    cbar.set_label(colorbarTitle, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)

    if toPlot:
        plt.show()

    if saveFig:
        fig.savefig(saveFileName)

    return fig, ax

def plotDiscreteBifurcationFile(inputFileName,
                saveFig = True, saveFileName="",
                transformData = lambda x:x,
                dicEqColor = {}, xlabel = "", ylabel="",
                tickText = [],
                toPlot=False, fontsize = 12):

    if saveFileName == "":
        saveFileName = inputFileName[0:-4] + ".png"

    data = pd.read_csv(inputFileName)
    data = data.to_numpy()
    listParamX = np.float64(data[0, 1:])
    listParamY = np.float64(data[1:, 0])
    color = np.array(data[1:, 1:], dtype=str)
    color = transformData(color)

    legendSize = 3
    myTickVals = [k for k in range(0,legendSize)]
    myTickBoundaries = [(2*k-1)/2 for k in range(0,legendSize+1)]

    if len(tickText) == 0:
        tickText = ["{0}".format(k) for k in range(legendSize)]

    fig, ax = plt.subplots(1)

    listColor = ["royalblue", "firebrick", "darkorange",
                    "green", "orchid", "black", "goldenrod",
                    "sienna", "purple", "navy", "silver", "turquoise"]
    counter = 0
    myColors = []

    for k in range(0, legendSize):
        if k in dicEqColor.keys():
            myColors.append(dicEqColor[k])
        else:
            myColors.append(listColor[counter])
            counter += 1
            print(k, " has no colour")

    print(counter, " colors were not defined")
    myColormap = ListedColormap(myColors)

    im = ax.pcolormesh(listParamX, listParamY, color, cmap=myColormap, vmin=myTickBoundaries[0], vmax=myTickBoundaries[-1])
    colorhatch = np.ma.masked_less(color, 2)

    ax.contourf(listParamX, listParamY, colorhatch, colors='none', hatches=['/'])
    cbar = fig.colorbar(im, boundaries = myTickBoundaries, values=myTickVals, drawedges = True)
    cbar.set_ticks(myTickVals)
    cbar.set_ticklabels(tickText, fontsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)

    if toPlot:
        plt.show()

    if saveFig:
        fig.savefig(saveFileName)

    return fig, ax

def plotSurface(inputFileNames,
                saveFig = True, saveFileName="",
                transformData = lambda x:x, xlabel = "", ylabel="",
                zlabel = "", colorbarLabels = ["", ""],
                toPlot=False, fontsize = 12):
    if saveFileName == "":
        saveFileName = inputFileNames[0][0:-4] + ".png"

    lambdaExistence = pd.read_csv(inputFileNames[0])
    lambdaExistence = lambdaExistence.to_numpy()
    listParamX = np.float64(lambdaExistence[0, 1:])
    listParamY = np.float64(lambdaExistence[1:, 0])
    valuesExistence = np.array(lambdaExistence[1:, 1:], dtype=float)

    lambdaStab = pd.read_csv(inputFileNames[1])
    lambdaStab = lambdaStab.to_numpy()
    valuesStab = np.array(lambdaStab[1:, 1:], dtype=float)


    fig, ax = plt.subplots(1, subplot_kw={"projection": "3d"})

    X, Y = np.meshgrid(listParamX, listParamY)

    myNorm = colors.LogNorm(vmin=min(valuesExistence.min(),valuesStab.min())
                            , vmax=max(valuesExistence.max(), valuesStab.max())
                            )
    myNorm1 = colors.LogNorm(vmin=valuesExistence.min()
                            , vmax=valuesExistence.max()
                            )
    # myNorm = colors.Normalize(vmin=min(valuesExistence.min(),valuesStab.min())
    #                         , vmax=max(valuesExistence.max(), valuesStab.max())
    #                         )


    colorStab = ax.plot_surface(X, Y, valuesStab, cmap=cm.autumn_r
                     ,norm = myNorm
                    )
    colorExistence = ax.plot_surface(X, Y, valuesExistence, cmap=cm.summer_r
                    ,norm = myNorm
                    )
    # cbarExistence = fig.colorbar(colorExistence, location='left', pad= 0.0)
    # cbarExistence.set_label(colorbarLabels[0], fontsize=fontsize)
    cbarStab = fig.colorbar(colorStab)
    cbarStab.set_label(colorbarLabels[1], fontsize=fontsize)
    ax.set_xlabel(xlabel,labelpad = 15, fontsize=fontsize)
    ax.set_ylabel(ylabel, labelpad = 15, fontsize=fontsize)
    ax.set_zlabel(zlabel, labelpad = 15, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)

    if toPlot:
        plt.show()

    if saveFig:
        fig.savefig(saveFileName)

    return fig, ax



def eqNamestoNbr(tab):
    rslt = np.zeros_like(tab, dtype=float)
    ncol, nrow = rslt.shape
    for indCol in range(ncol):
        for indRow in range(nrow):
            val = tab[indCol, indRow]
            if val == "HFW":
                rsltVal = 1
            elif val == "LC":
                rsltVal = 2
            else: ## NaN case
                rsltVal = 0
            rslt[indCol, indRow] = rsltVal
    return rslt

def addHDHW(tab):
    nrow, _ = tab.shape
    rslt = np.zeros((nrow, 2), dtype=float)
    for indRow in range(nrow):
            H = tab[indRow, 1] + tab[indRow, 3]
            FW = tab[indRow, 2]
            rslt[indRow, 0] = H
            rslt[indRow, 1] = FW
            
    return rslt


input = ["/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/GradI/inputI0/result.csv"
        ,"/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/GradI/inputI1/result.csv"
        ,"/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/GradI/inputI10/result.csv"
        # ,"/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/FromEEFtoEEHF/inputI5/result.csv"
        # ,"/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/FromEEFtoEEHF/inputI7_5/result.csv"
        # ,"/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/FromEEFtoEEHF/inputI10/result.csv",
        # "/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/FromEEFtoEEHF/input7/result.csv"
        ]


myFig, myAx = plotPhasePlane(input, toPlot=False,
                                transformData=addHDHW,
                                xlabel = "$H_D + H_W$", ylabel="$F_W$",
                               fontsize=20, markersize = 40, lineColors=['red', 'blue', 'black'],
                               lineLabels=["$0$", "1", "5.5", "$5$", "7.5", "$10$", "150"],
                               titleLegend="$\\mathcal{I} = $")

plt.show(block = True)




