import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib import ticker

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
                lineLabels = [], titleLegend = ""):
    
    if saveFileName == "":
        saveFileName = inputFileNames[0][0:-4] + ".png"
    if len(lineLabels) < len(inputFileNames):
        lineLabels = lineLabels + ["" for _ in range(len(inputFileNames)- len(lineLabels))]

    fig, ax = plt.subplots(1)

    for k, inputFileName in enumerate(inputFileNames):

        data = pd.read_csv(inputFileName)
        data = data.to_numpy()
        dataTransformed = transformData(data)
        print(data.shape)
        listParamX = np.float64(dataTransformed[1:, 0])
        listParamY = np.float64(dataTransformed[1:, 1])
        line, = ax.plot(listParamX, listParamY, label=lineLabels[k])
        ax.scatter(listParamX[0], listParamY[0], marker='+', s=markersize, c=line.get_color())
        ax.scatter(listParamX[-1], listParamY[-1], marker='o', s=markersize, c=line.get_color())


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

def plotSurfaceTest(inputFileNames,
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


    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, subplot_kw={"projection": "3d"}, sharex=True)

    X, Y = np.meshgrid(listParamX, listParamY)

    myNorm = colors.LogNorm(vmin=min(valuesExistence.min(),valuesStab.min())
                            , vmax=max(valuesExistence.max(), valuesStab.max())
                            )
    myNormExistence = colors.LogNorm(vmin=valuesExistence.min()
                            , vmax=valuesExistence.max()
                            )
    myNormStab = colors.LogNorm(vmin=valuesStab.min()
                            , vmax=valuesStab.max()
                            )

    colorExistence = ax1.plot_surface(X, Y, valuesExistence, cmap=cm.summer_r
                    ,norm = myNormExistence
                    )
    colorExistence2 = ax2.plot_surface(X, Y, valuesExistence, cmap=cm.summer_r
                    ,norm = myNormExistence
                    )
    colorStab = ax2.plot_surface(X, Y, valuesStab, cmap=cm.autumn_r
                     ,norm = myNormStab
                    )
    cbarExistence = fig.colorbar(colorExistence, location='left', pad= 0.0)
    cbarExistence.set_label(colorbarLabels[0], fontsize=fontsize)
    cbarStab = fig.colorbar(colorStab)
    cbarStab.set_label(colorbarLabels[1], fontsize=fontsize)

    ax1.set_xlabel(xlabel,labelpad = 15, fontsize=fontsize)
    ax1.set_ylabel(ylabel, labelpad = 15, fontsize=fontsize)
    ax1.set_zlabel(zlabel, labelpad = 15, fontsize=fontsize)
    ax1.tick_params(labelsize=fontsize)

    ax2.set_xlabel(xlabel,labelpad = 15, fontsize=fontsize)
    ax2.set_ylabel(ylabel, labelpad = 15, fontsize=fontsize)
    ax2.set_zlabel(zlabel, labelpad = 15, fontsize=fontsize)
    ax2.tick_params(labelsize=fontsize)


    if toPlot:
        plt.show()

    if saveFig:
        fig.savefig(saveFileName)

    return fig, (ax1, ax2)

def myTransformData(tab):
    rslt = np.zeros_like(tab, dtype=float)
    ncol, nrow = rslt.shape
    for indCol in range(ncol):
        for indRow in range(nrow):
            val = tab[indCol, indRow]
            if val > 0:
                rsltVal = 2
            elif val < 0:
                rsltVal = 1
            else: ## NaN case
                rsltVal = 0
            rslt[indCol, indRow] = rsltVal
    return rslt

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

# input = "/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/stabcI0/Plot3/BifurcationcIM02/Bifurcation.csv"

# myFig, myAx = plotDiscreteBifurcationFile(
#                     input, ylabel="$K_F(1-\\alpha)$", xlabel="$\\lambda_{F}$",
#                     saveFig = False, toPlot=False, fontsize=20
#                     , tickText=["$EE^{F}$ AS", "$EE^{HF_W}$ AS", "$EE^{HF_W}$ \n Limit Cycle"]
#                     , dicEqColor={0 : "white", 1 : "grey", 2 : "black"}
#                     , transformData=eqNamestoNbr
#                     )

# inputLambdaMax = "/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/stabcI0/Plot3/BifurcationcIM02/LambdaMax.csv"
# data = pd.read_csv(inputLambdaMax, low_memory=False)
# data = data.to_numpy()
# listAlpha = np.float64(data[1:, 0])
# listLambdaMax = np.float64(data[1:, 1])
# line, = myAx.plot(listLambdaMax, listAlpha, '--', linewidth = 3,
#            color = "orange")
# line.set_label("$\\lambda_{F, I=0}^{Max}(K_{F,\\alpha})$")

# inputLambdaMin = "/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/stabcI0/Plot3/BifurcationcIM02/LambdaMin.csv"
# data = pd.read_csv(inputLambdaMin, low_memory=False)
# data = data.to_numpy()
# listAlpha = np.float64(data[1:, 0])
# listLambdaMin = np.float64(data[1:, 1])
# line, = myAx.plot(listLambdaMin, listAlpha, '--', linewidth = 3,
#            color = "green")
# line.set_label("$\\lambda_{F, I=0}^{Min}(K_{F,\\alpha})$")
# myAx.legend(fontsize=20)
# myAx.set_xlim((0.0, 0.2))

# line, = myAx.plot(2.304 * np.ones(50), np.linspace(min(listAlpha), max(listAlpha), 50),
#                   '--', linewidth = 3,
#                     color = "red")
# line.set_label("$\\lambda_{F, I=0.1}^{Max}(K_{F,\\alpha})$")
# myAx.legend(fontsize=20)

# # myAx.set_xlim((0.0, 0.75))

# # line, = myAx.plot(0.1918* np.ones(100), 3000*(1- np.linspace(0, 1, 100)), '--', color = "red",
# #                     linewidth=3)
# # line.set_label("$\\lambda_F^{Max}$")
# # myAx.set_ylim(bottom=0)
# # myAx.legend(fontsize=20)

# input3D = ["/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/stabcI0/Plot3/lambdaMAlphaExistence.csv"
#            ,"/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/stabcI0/Plot3/lambdaMAlphaStab.csv"
#             ]

# Fig, myAx = plotSurface(input3D
#                       ,toPlot=False
#                       ,saveFig= False
#                       ,ylabel="$K_F(1-\\alpha)$"
#                       ,xlabel="$m$"
#                       ,zlabel = "$\\lambda_F$"
#                       ,colorbarLabels=["$\\lambda_{F, I=0}^{Min}(m; K_{F,\\alpha})$", "$\\lambda_{F, I=0}^{Max}(m; K_{F,\\alpha})$"]
#                       ,fontsize=20
#                       )
# myAx.locator_params(axis='y', nbins = 3)
# myAx.locator_params(axis='z', nbins = 3)


input = ["/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/FromEEFtoLC/EEF/result.csv",
                  "/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/FromEEFtoLC/EELC/result.csv",
         "/home/hetier/Documents/0PhD/SimulationNumerique/HunterRC/FromEEFtoLC/EEHF/result.csv",

        ]


myFig, myAx = plotTrajectory1D(input, toPlot=False,
                            #    xlabel="$H_D$", ylabel="$F_W$", zlabel="$H_W$",
                                xlabel="time (year)", y1label="$H_D$", y2label="$F_W$", y3label="$H_W$",
                               fontsize=20, markersize = 40,
                               lineLabels=["$0$", "1", "2.5", "$5$", "7.5", "$10$", "150"],
                               titleLegend="$\\mathcal{I} = $")

# myAx.locator_params(axis='y', nbins = 5)
myAx3 = myAx[2]
myAx2 =myAx[1]
myAx1 = myAx[0]
# myAx1.set_xlim((0,60))

box = myAx1.get_position()
myAx1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
box = myAx2.get_position()
myAx2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
box = myAx3.get_position()
myAx3.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
myAx2.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20, 
                title_fontsize = 20, title = "$\\mathcal{I} = $")


plt.show(block = True)




