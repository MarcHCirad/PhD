import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib import ticker
import pathlib


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

    legendSize = len(tickText)
    myTickVals = [k for k in range(0,legendSize)]
    myTickBoundaries = [(2*k-1)/2 for k in range(0,legendSize+1)]

    if len(tickText) == 0:
        tickText = ["{0}".format(k) for k in range(legendSize)]

    fig, ax = plt.subplots(1)

    listColor = ["royalblue", "firebrick", "darkorange"]
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
    # colorhatch = np.ma.masked_less(color, 3)

    ax.contourf(listParamX, listParamY, color, colors='none')
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
            elif val == "H":## NaN case
                rsltVal = 0
            elif val == "F":## NaN case
                rsltVal = 0
            elif val == "HHFW2?":
                rsltVal = 3
            else:
                rsltVal = 4

            rslt[indCol, indRow] = rsltVal
    return rslt

pathWD = str(pathlib.Path(__file__).parent.resolve())

input = pathWD + "/Diagram.csv"

color0 = "#ea801c"
color0 = "#A00000"
color2 = "#C9FFFF"
fontsize = 30

myFig, myAx = plotDiscreteBifurcationFile(
                    input, ylabel="$\\alpha$", xlabel="$\\lambda_{F}$",
                    saveFig = False, toPlot=False, fontsize=fontsize
                    , tickText=["$EE^{F}$ GAS", "$EE^{HF_W}$ GAS", "Limit Cycle  \n  around $EE^{HF_W}$"]
                    , dicEqColor={0 : "#A00000", 1 : "#298c8c", 2 : "#C9FFFF"}#, 3:"green", 4 :"red"}
                    , transformData=eqNamestoNbr
                    )

xlim = myAx.get_xlim()

inputLambdaMin = pathWD + "/Diagram1DMin.csv"
data = pd.read_csv(inputLambdaMin, low_memory=False)
data = data.to_numpy()
listAlpha = np.float64(data[1:, 0])
listLambdaMin = np.float64(data[1:, 1])
line, = myAx.plot(listLambdaMin[0:-1:10], listAlpha[0:-1:10], '--', linewidth = 5,
           color = "black")
line.set_label("$\\lambda_{F, I=0}^{Min}(\\alpha)$")
myAx.legend(fontsize=fontsize)


inputLambdaMax = pathWD + "/Diagram1DMax.csv"
data = pd.read_csv(inputLambdaMax, low_memory=False)
data = data.to_numpy()
listAlpha = np.float64(data[1:, 0])
listLambdaMax = np.float64(data[1:, 1])
line, = myAx.plot(listLambdaMax[0:-1:10], listAlpha[0:-1:10], '-.', linewidth = 3,
           color = "black")
line.set_label("$\\lambda_{F, \mathcal{I}=0}^{Max}(\\alpha)$")
myAx.legend(fontsize=fontsize, loc="upper right", title="$\mathcal{I} = 0$", title_fontsize = fontsize)
myAx.set_xlim(xlim)
# myAx.set_xlim((0, 1e-5))

# myAx.xaxis.get_major_formatter().set_scientific(False)
# myAx.set_xticklabels(['0', '2e-6', '4e-6', '6e-6','8e-6','1e-5'])


plt.show(block = True)



