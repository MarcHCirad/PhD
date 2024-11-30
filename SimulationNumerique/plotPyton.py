import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors

def plotPhasePortrait(inputFileNames, variables,
                    saveFig = True,
                    saveFileName = "",
                    title = "",
                    legend =  {},
                    listColor = [],
                    toPlot = False,
                    fontsize = 12,
                    transform = ""):
    inputFileNames = [name + "/result.csv" for name in inputFileNames]
    if not listColor:
        listColor = ["royalblue", "firebrick", "darkorange", 
                     "green", "orchid", "black", "goldenrod",
                    "sienna", "purple", "navy", "silver", "turquoise"]
    
    if saveFileName == "":
        saveFileName = inputFileNames[1][0:-4] + ".png"

    fileName = inputFileNames[0]
    with open(fileName, newline='') as myCSV:
        spamReader = csv.reader(myCSV, delimiter=",")
        header = next(spamReader)

        indVariable1 = variables[0] + 1
        indVariable2 = variables[1] + 1

        if transform == "plus":
            nameVariable2 = "F_D + V_D"
        else:
            nameVariable2 = header[indVariable2]
        
        nameVariable1 = header[indVariable1]
        
        dataVariable1, dataVariable2 = [], []
        for row in spamReader:
            dataVariable1.append(float(row[indVariable1]))
            
            if transform == "plus":
                dataVariable2.append(float(row[indVariable2]) + float(row[2]))
            else:
                dataVariable2.append(float(row[indVariable2]))
    if 1 in list(legend):
        showLegend = True
        legendText = legend[1]
    else:
        showLegend = False
        legendText = ""
    
    myFig, myAxes = plt.subplots(1)
    line, = myAxes.plot(dataVariable1, dataVariable2, marker="", color=listColor[0])
    if showLegend:
        line.set_label(legendText, fontsize=fontsize)

    myAxes.plot(dataVariable1[0], dataVariable2[0], '+', color=listColor[0], markeredgewidth= fontsize/4)
    myAxes.plot(dataVariable1[-1], dataVariable2[-1], 'o', color=listColor[0], markeredgewidth= fontsize/2)

    for (k, fileName) in enumerate(inputFileNames[1:]):
        with open(fileName, newline='') as myCSV:
            spamReader = csv.reader(myCSV, delimiter=",")
            _ = next(spamReader)           
            dataVariable1, dataVariable2 = [], []
            for row in spamReader:
                dataVariable1.append(float(row[indVariable1]))
                if transform == "plus":
                    dataVariable2.append(float(row[indVariable2]) + float(row[2]))
                else:
                    dataVariable2.append(float(row[indVariable2]))

        if k+2 in list(legend):
            showLegend = True
            legendText = legend[k+2]
        else:
            showLegend = False
            legendText = ""
        line, = myAxes.plot(dataVariable1, dataVariable2, marker="", color=listColor[k+1])
        
        if showLegend:
            line.set_label(legendText, fontsize = fontsize)
            
        myAxes.plot(dataVariable1[0], dataVariable2[0], '+', color=listColor[k+1], markeredgewidth= fontsize/4)
        myAxes.plot(dataVariable1[-1], dataVariable2[-1], 'o', color=listColor[k+1], markeredgewidth= fontsize/2)


    myAxes.set_title(title)
    if showLegend:
        myAxes.legend()
    myAxes.set_xlabel("$" + nameVariable1 +"$", fontsize = fontsize)
    myAxes.set_ylabel("$" + nameVariable2 +"$", fontsize = fontsize)
    myAxes.tick_params(axis='both', which='major', labelsize=fontsize)

    if toPlot:
        myFig.show()
        plt.waitforbuttonpress()

    if saveFig:
        myFig.savefig(saveFileName)
    
    return myFig, myAxes


def transformEqNames(listNames):
    rslt = []
    for name in listNames:
        if name == "(TE)":
            newName = "TE"
        else:
            newName = name.replace("(", "EE^{")
            newName = newName.replace(")", "}, ")
            newName = newName.replace("A", "D")
            newName = newName.replace(",3}", "}_3")
            newName = newName.replace("EE^{TE}", "TE")
            newName = newName[:-2]

        newName = "$" + newName + "$"
        rslt.append(newName)
    return rslt

def plotBifurcationFile(inputFileName,
                        saveFig = True, saveFileName="",
                        dicEqColor = {}, xlabel = "", ylabel="",
                        toPlot=False, fontsize = 12):
    
    if saveFileName == "":
        saveFileName = inputFileName[0:-4] + ".png"

    with open(inputFileName, newline='') as myCSV:
        spamReader = csv.reader(myCSV, delimiter=",")
        _ = next(spamReader)
        myRow = next(spamReader)
        listParamY = [float(val) for val in myRow[1:]]
        
        dimY = len(listParamY)
        dimX = sum(1 for _ in csv.reader(myCSV))

    with open(inputFileName, newline='') as myCSV:
        spamReader = csv.reader(myCSV, delimiter=",")
        _ = next(spamReader)
        _ = next(spamReader)
        listParamX = []
               
        eqNames = np.zeros((dimX, dimY), dtype=object)
        for (k, row) in enumerate(spamReader):
            listParamX.append(float(row[0]))
            eqNames[k, : ] = row[1:]


    eqNamesUnique = np.unique(eqNames)
    legendSize = len(eqNamesUnique)


    dicEqNbr = {eqNamesUnique[i]:i for i in range(0,legendSize)}
    dicNbrEq = {i:eqNamesUnique[i] for i in range(0,legendSize)}

    color = np.zeros((dimY, dimX))

    for ind_col in range(dimY):
        for ind_row in range(dimX):
            color[ind_col, ind_row] = dicEqNbr[eqNames[ind_row, ind_col]]
 
    

    myTickVals = [k for k in range(0,legendSize)]
    myTickBoundaries = [(2*k-1)/2 for k in range(0,legendSize+1)]

    # myTickText = np.sort(dicEqNbr.keys())
    myTickText = transformEqNames(dicEqNbr.keys())
    
    fig, ax = plt.subplots(1)

    listColor = ["teal", "cyan", "red", "black", "ForestGreen"]
    counter = 0
    myColors = []
    for k in range(0, legendSize):
        if dicNbrEq[k] in dicEqColor.keys():
            myColors.append(dicEqColor[dicNbrEq[k]])
        else:
            myColors.append(listColor[counter])
            counter += 1
            print(dicNbrEq[k], " has no colour")
    
    print(counter, " colors were not defined")

    myColormap = ListedColormap(myColors)

    im = ax.pcolormesh(listParamX, listParamY, color, cmap=myColormap)
    cbar = fig.colorbar(im, boundaries = myTickBoundaries, values=myTickVals) 
    cbar.set_ticks(myTickVals)
    cbar.set_ticklabels(myTickText, fontsize=fontsize)
    ax.set_xlabel("$" + xlabel + "$", fontsize=fontsize)
    ax.set_ylabel("$" + ylabel + "$", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    
    if toPlot:
        plt.show()
    
    if saveFig:
        fig.savefig(saveFileName)
    
    return fig, ax

def plotBifurcationFileV2(inputFileName,
                        saveFig = True, saveFileName="",
                        dicEqColor = {}, xlabel = "", ylabel="",
                        toPlot=False, fontsize = 12):
    
    if saveFileName == "":
        saveFileName = inputFileName[0:-4] + ".png"

    data = pd.read_csv(inputFileName)
    data = data.to_numpy()
    listParamX = np.array(data[0, 1:], dtype=float)
    listParamY = np.array(data[1:, 0], dtype=float)

    eqNames = np.array(data[1:, 1:], dtype=float)

    eqNamesUnique = np.unique(eqNames)
    legendSize = len(eqNamesUnique)


    dicEqNbr = {eqNamesUnique[i]:i for i in range(0,legendSize)}
    dicNbrEq = {i:eqNamesUnique[i] for i in range(0,legendSize)}

    color = np.zeros_like(eqNames)
    nbrRow, nbrCol = color.shape

    for ind_col in range(nbrCol):
        for ind_row in range(nbrRow):
            color[ind_row, ind_col] = dicEqNbr[eqNames[ind_row, ind_col]]
 
    myTickVals = [k for k in range(0,legendSize)]
    myTickBoundaries = [(2*k-1)/2 for k in range(0,legendSize+1)]
    myTickText = transformEqNames(dicEqNbr.keys())
    
    fig, ax = plt.subplots(1)

    listColor = ["teal", "cyan", "red", "black", "ForestGreen"]
    counter = 0
    myColors = []
    for k in range(0, legendSize):
        if dicNbrEq[k] in dicEqColor.keys():
            myColors.append(dicEqColor[dicNbrEq[k]])
        else:
            myColors.append(listColor[counter])
            counter += 1
            print(dicNbrEq[k], " has no colour")
    
    print(counter, " colors were not defined")

    myColormap = ListedColormap(myColors)

    im = ax.pcolormesh(listParamX, listParamY, color, cmap=myColormap)
    cbar = fig.colorbar(im, boundaries = myTickBoundaries, values=myTickVals) 
    cbar.set_ticks(myTickVals)
    cbar.set_ticklabels(myTickText, fontsize=fontsize)
    ax.set_xlabel("$" + xlabel + "$", fontsize=fontsize)
    ax.set_ylabel("$" + ylabel + "$", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    
    if toPlot:
        plt.show()
    
    if saveFig:
        fig.savefig(saveFileName)
    
    return fig, ax


def plotValEqFile(inputFileName,
                saveFig = True, saveFileName="",
                xlabel = "", ylabel="", colorbarTitle = "",
                toPlot=False, fontsize = 12):
    
    if saveFileName == "":
        saveFileName = inputFileName[0:-4] + ".png"

    data = pd.read_csv(inputFileName)
    data = data.to_numpy()
    listParamX = np.float64(data[0, 1:])
    listParamY = np.float64(data[1:, 0])
    values = np.array(data[1:, 1:], dtype=float)
                     
    fig, ax = plt.subplots(1)

    im = ax.pcolormesh(listParamX, listParamY, values
                       ,norm=colors.LogNorm(vmin=values.min(), vmax=values.max())
                    #    ,norm = colors.PowerNorm(gamma=0.05)
                    #    ,norm=colors.SymLogNorm(linthresh = 0.0001, vmin=values.min(), vmax=values.max())
                    )
    cbar = fig.colorbar(im)
    cbar.set_label(colorbarTitle, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    
    ax.set_xlabel("$" + xlabel + "$", fontsize=fontsize)
    ax.set_ylabel("$" + ylabel + "$", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    
    if toPlot:
        plt.show()
    
    if saveFig:
        fig.savefig(saveFileName)
    
    return fig, ax


def plotStabFile(inputFileName,
                saveFig = True, saveFileName="",
                dicEqColor = {}, xlabel = "", ylabel="",
                toPlot=False, fontsize = 12):
    
    if saveFileName == "":
        saveFileName = inputFileName[0:-4] + ".png"

    data = pd.read_csv(inputFileName)
    data = data.to_numpy()
    listParamX = np.float64(data[0, 1:])
    listParamY = np.float64(data[1:, 0])
    color = np.array(data[1:, 1:], dtype=float)
               
    legendSize = 3
    myTickVals = [k for k in range(0,legendSize)]
    myTickBoundaries = [(2*k-1)/2 for k in range(0,legendSize+1)]
    myTickText = ["$\\Delta_{Stab}(\\alpha, \\lambda_{F}) < 0$", "$\\Delta_{Stab}(\\alpha, \\lambda_{F}) = 0$",  "$\\Delta_{Stab}(\\alpha, \\lambda_{F}) > 0$"]
    
    fig, ax = plt.subplots(1)

    listColor = ["teal", "cyan", "red", "black", "ForestGreen"]
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

    im = ax.pcolormesh(listParamX, listParamY, color, cmap=myColormap)
    cbar = fig.colorbar(im, boundaries = myTickBoundaries, values=myTickVals)
    cbar.set_ticks(myTickVals)
    cbar.set_ticklabels(myTickText, fontsize=fontsize)
    ax.set_xlabel("$" + xlabel + "$", fontsize=fontsize)
    ax.set_ylabel("$" + ylabel + "$", fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    
    if toPlot:
        plt.show()
    
    if saveFig:
        fig.savefig(saveFileName)
    
    return fig, ax



# input = "/home/hetier/Documents/PhD/SimulationNumerique/HunterTest/stab/stabDiagramLambdaRealTest.csv"
# myFig, myAx = plotValEqFile(input, ylabel="\\alpha", xlabel="\\lambda_{F}", saveFig = False,
#                     toPlot=False, fontsize=20
#                     ,colorbarTitle="$SymLogNorm(\\Delta_{stab})$"
#                     # , dicEqColor={0 : "white", 1 : "grey", 2 : "black"}
#                     )

# lim = myAx.get_xlim()
# myAx.set_xticks(list(myAx.get_xticks()) + [0.1919])
# myAx.set_xlim(lim)
# tickname = [name for name in myAx.get_xticklabels()]
# tickname[-1] = "$\\lambda_{F, max}$"
# myAx.set_xticklabels(tickname)
# myFig.show()
# myFig.waitforbuttonpress()

# inputFeq = "/home/hetier/Documents/PhD/SimulationNumerique/HunterTest/stab/FeqDiagramLambdaRealTest.csv"
# myFig, myAx = plotValEqFile(inputFeq, ylabel="\\alpha", xlabel="\\lambda_{F}", saveFig = False,
#                     toPlot=False, fontsize=20
#                     ,colorbarTitle="$F^*_W$"
#                     )
# lim = myAx.get_xlim()
# myAx.set_xticks(list(myAx.get_xticks()) + [0.1919])
# myAx.set_xlim(lim)
# tickname = [name for name in myAx.get_xticklabels()]
# tickname[-1] = "$\\lambda_{F, max}$"
# myAx.set_xticklabels(tickname)
# myFig.show()
# myFig.waitforbuttonpress()

# input = "/home/hetier/Documents/PhD/SimulationNumerique/Hunter/bifurcation/bifurcationDiagramCLambdaFWH.csv"

# dicEqNamesColor = {"(F_W)" : "silver",
#                     "(F_W)(H_AH_W)" : "firebrick",
#                     "(F_W)(H_AF_WH_W)" : "green"}

# myFig, myAx = plotBifurcationFile(input, xlabel="\\lambda_{F_D}", ylabel="\\lambda_{V_D}", saveFig = False,
#                     toPlot=False, fontsize=20, dicEqColor=dicEqNamesColor)



# line, = myAx.plot(np.linspace(2.9, 10), 4*0.64/(200 * 0.2/0.5 * 1)* np.ones(50), '--', linewidth = 5,
#            color = "orange")
# line.set_label("$T_F(a_W K_F) = 1/2$")

# TFc = 0.64 / (0.2/0.5) / np.linspace(2.9, 10, 200)

# line2, = myAx.plot(np.linspace(2.9, 10, 200), TFc, color = "blue", linewidth = 5)
# line2.set_label("$T_F(c) = 1$")

# myAx.legend(fontsize=20)

# plt.show(block = True)
# myFig.savefig(input[0:-4] + "AddTFc.png")


inputNames = ["/home/hetier/Documents/PhD/SimulationNumerique/HunterTest/hunterRK4_"
            #   ,"/home/hetier/Documents/PhD/SimulationNumerique/HunterTest/hunterRK4_1"
              ]

saveFileName="/home/hetier/Documents/PhD/SimulationNumerique/HunterTest/Orbit.png"

myFig, myAx = plotPhasePortrait(inputNames, [0, 1], listColor=["red",  "blue"],
                                transform="plusa", fontsize=20, saveFileName=saveFileName)

myAx.plot([12.48], [159.83], marker = '+', linestyle = "--", color="black", markeredgewidth = 4)
# myAx.set_ybound(lower=-4, upper=90)
# lim = myAx.get_xlim()
# myAx.set_xticks(list(myAx.get_xticks()) + [25])
# myAx.set_xlim(lim)
# tickname = [name for name in myAx.get_xticklabels()]
# tickname[6] = "$\\beta$"
# myAx.set_xticklabels(tickname)

myFig.show()
myFig.waitforbuttonpress()
myFig.savefig(saveFileName[0:-4] + "TT.png")


