t = @elapsed begin

abstract type mathematicalModel end
abstract type numericalModel end

include("modelWild.jl")
include("modelWildRK4.jl")
include("modelAnthropized.jl")
include("modelAnthropizedRK4.jl")
include("modelMigration.jl")
include("modelMigrationRK4.jl")

include("plotFile.jl")
include("writer.jl")
include("reader.jl")



function solveModel(allModels::Dict{String, numericalModel})
    for model in collect(values(allModels)) 
        solveModel(model)
    end
end

function solveWriteModel(inputFile::String)
    myModel = readNumericalModel(inputFile)
    solveModel(myModel)
    indInput = first(findfirst("input", inputFile))
    dirName = inputFile[1:last(indInput)-2] ## we delete .txt
    suffix = inputFile[first(indInput)+5:end-4]
    resultFolder = writeNumericalModel(myModel, dirName, suffix)
    cp(inputFile, resultFolder * "/" * inputFile[first(indInput):end], force=true)
    return resultFolder
end

function solveWriteModels(dirName::String)
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    inputNames = searchdir(dirName, "input")

    resultFolders = []
    for inputName in inputNames
        push!(resultFolders, solveWriteModel(dirName * inputName))
    end
    return resultFolders
end

function writeNumericalModel(allModels::Dict{String, numericalModel}, dirPrefix::String)
    listFolderNames = Vector{String}()
    for model in collect(values(allModels))
        push!(listFolderNames, writeNumericalModel(model, dirPrefix))
    end
    return listFolderNames
end

function readCreateNumericalModel(inputFile::String, nameParam::String, listParam::Vector{Float64})
    fileNames = []
    for (ind, param) in enumerate(listParam)
        name = inputFile[1:end-4] * nameParam * string(ind) * ".txt"
        append!(fileNames, name)
        open(inputFile, "r") do fin
            open(name, "w") do fout
            for line in eachline(fin)
                if contains(line, nameParam)
                    write(fout, nameParam * "    " * string(param) *"\n")
                else
                    write(fout, line *"\n")
                end
            end
            end
        end
    end
    return fileNames
end

function mainBifurcation()

    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"
    pathAnthropo = pathWD * "/Anthropized/"

    myModel = readNumericalModel(pathAnthropo * "/input.txt")
    math = myModel.mathModel  
    
    ## color ##
    dicEqNamesColor = Dict(["(F)" => "powderblue",
                        "(F), (FH)" => "firebrick",
                        "(F), (FH), (FVH_2)" => "turquoise",
                        "(F), (FH), (VH_2)" => "darkorange",
                        "(F), (FH), (VH_2), (FVH_2)" => "gold",
                        "(F), (FH_\beta)" => "green",
                        "(F), (FH_\beta), (FH)" => "khaki",
                        "(F), (FH_\beta), (FH), (FVH_2)" => "skyblue",
                        "(F), (FH_\beta), (FH), (VH_2)" => "tan",
                        "(F), (FH_\beta), (FH), (VH_2), (FVH_2)" => "chocolate",
                        "(F), (FH_\beta), (VH_2)" => "olivedrab",
                        "(F), (FVH_\beta)" => "silver",
                        "(F), (FVH_\beta), (FH)" => "indianred",
                        "(F), (FVH_\beta), (FH), (FVH_2)" => "black",
                        "(F), (FVH_\beta), (FH), (VH_2)" => "navy",
                        "(F), (FVH_\beta), (FH), (VH_2), (FVH_2)" => "sienna",
                        "(F), (FVH_\beta), (FVH_2)" => "royalblue",
                        "(F), (FVH_\beta), (VH_2)" => "goldenrod",
                        "(F), (FVH_\beta), (VH_2), (FVH_2)" => "seashell",
                        "(F), (H)" => "gold",
                        "(F), (H), (VH_2)" => "bisque",
                        "(F), (H_\beta)" => "darkgreen",
                        "(F), (VH_\beta)" => "seagreen",
                        "(F), (VH_\beta), (VH_2)" => "purple",
                        "(F), (VH_2)" => "lightgreen"
                    ])
    dicEqNamesColor = Dict(["(TE)(H_A)" => "silver",
                            "(TE)(H_A)(V_AH_A)" => "firebrick",
                            "(TE)(V_AH_A)" => "green",
                            "(V_A)(H_A)" => "seashell",
                            "(V_A)(V_AH_A)" => "sienna",
                            "(TE)(V_AH_A_1)(V_AH_A)" => "skyblue"])

    
    
    ## parameter list ##
    listLambdaFH = [k for k = 0.001:0.01:0.5]
    listLambdaVH = [k for k = 0.001:0.0001:0.05]
    listLV = [k for k = 1.0:0.1:75.1]
    listbeta = [k for k = 0.1:0.01:12]
    lista = [k for k = 0.01:0.05:15]
    listb = [k for k = 0.01:0.05:15]
    listc = [k for k = 2.9:0.1:3.1]

    ## dir and file name ##
    
    dirName = pathAnthropo * "/bifurcation/comparisonV2/"
    nameSuffix = "bifurcationDiagramCLambdaVH3"
    println(math)

    ## computation ##
    bifurcationMatrix = computeBifurcationDiagramFH(math, "c", listc, "lambdaFH", listLambdaVH ; eqValues=false)
    filePrefix = writeBifurcationDiagram(math, bifurcationMatrix, dirName; 
                                        nameSuffix = nameSuffix,
                                        eqVals = false)
    ## plot ##
    plotBifurcationFile(filePrefix * ".csv", 
                        filePrefix * ".html" ;
                        # title = latexstring("\\text{Bifurcation Diagram in } (c, L_V) \\text{ plane}"),
                        title = "",
                        xlabel = latexstring("\\LARGE{c}"),
                        ylabel = latexstring("\\LARGE{\\lambda_{VH}}"),
                        dicEqColor = dicEqNamesColor,
                        toPlot = true,
                        eqVals = false,
                        titleEqVals = "\\text{: average of equilibrium values}",
                        fontsize = 30)
end

function mainSolveModel()
    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"
    pathWild = pathWD * "/Wild/"
    pathAnthropo = pathWD * "/Anthropized/"
    pathMigration = pathWD * "/Migration/"
    path = pathAnthropo
    
    # resultFolders = ["/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_",
                    # "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_2", 
                    # "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_3", 
                    # "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_4", 
                    # "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_5",
                    # "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_6",
                    # "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_7"]
    
    # resultFolders = solveWriteModels(path)
    # myModel = readNumericalModel(path * "/input.txt")
    # println("Eq FV : ", equilibriumFV(myModel.mathModel))
    # println("Eq HFV : ", equilibriumHFV(myModel.mathModel))
    # println(computeThresholdsFH(myModel.mathModel))
    # println(interpretThresholdsFH(myModel.mathModel, computeThresholdsFH(myModel.mathModel)))


    ## Plotting ##
    # println([resultFolder * "/result.csv" for resultFolder in resultFolders])

    resultFolders = ["/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_01/", 
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_02/", 
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_04/", 
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_05/", 
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_06/", 
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_07/", 
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_10/", 
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_11/",
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_12/",
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_13/", 
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_21/",
    "/home/hetier/Documents/PhD/SimulationNumerique/Anthropized/anthropizedRK4_22/"]

    listColor = ["red", "red", "red", "red", "red", "red", "blue", "blue", "blue", "blue",
                "green", "green", "green"]
    plotPhasePortrait([resultFolder * "/result.csv" for resultFolder in resultFolders],
                        path*"/plotOrbit.html",
                        [1,3];
                        # title=latexstring("\\text{Orbit in the } F_W - V_W \\text{ plane}"),
                        legend=Dict([1=>latexstring("EE^{F_AH_A}_1"), 7=>latexstring("EE^{F_AH_A}_1"), 11=>latexstring("TE")]),
                        listColor = listColor,
                        completeLegend = true,
                        showLegend = true,
                        toPlot=true,
                        fontsize = 30)
    # plotPhasePortrait([resultFolder * "/result.csv" for resultFolder in resultFolders],
    #                     pathMigration*"/plotHF.html",
    #                     [1,2];
    #                     # title=latexstring("\\text{Orbit in the } F_W - V_W \\text{ plane}"),
    #                     # legend=[latexstring("allee"), "ecoService"],
    #                     showLegend = false,
    #                     toPlot=true,
    #                     fontsize = 30)
    # plotTrajectory1d([resultFolder * "/result.csv" for resultFolder in resultFolders],
    #                     path*"/plotHFV1D.html",
    #                     3;
    #                     # title=latexstring("\\text{Orbit in the } F_W - V_W \\text{ plane}"),
    #                     # legend=[latexstring("allee"), "ecoService"],
    #                     showLegend = true,
    #                     toPlot=true,
    #                     fontsize = 30)
end

function comparisonH()
    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"

    pathWD = pathWD * "/TestA"
    allModels = readNumericalModel(pathWD)
    numModel = allModels["Anthropized"]
    mathModel = numModel.mathModel

    listLV = [k for k = 0.1:0.1:70]
    listLF = [k for k = 0.1:0.1:70]
    mylist = [k for k = 0.1:0.01:7]

    ## dir and file name ##
    dirName = pathWD * "/comparisonH/"
    nameSuffix = "bifurcationDiagramCLV"
    filePrefix = dirName * "/" * nameSuffix

    valFH = []
    valFVH = []

    localParam = Dict([
        ("rF",mathModel.rF), ("KF",mathModel.KF), ("LF",mathModel.LF), ("muF",mathModel.muF), ("lambdaFH", mathModel.lambdaFH),
        ("rV",mathModel.rV), ("KV",mathModel.KV), ("LV", mathModel.LV), ("muV",mathModel.muV), ("lambdaVH", mathModel.lambdaVH),
        ("rH",mathModel.rH), ("a", mathModel.a), ("b", mathModel.b), ("c", mathModel.c), ("beta", mathModel.beta)
        ])

    mylist = listLF
    paramName = "LF"
    
    for param in mylist
        localParam[paramName] = param
        localModel = modelAnthropized(localParam)
        append!(valFH, equilibriumFH(localModel)[3])
        append!(valFVH, equilibriumFVH(localModel)[3])
    end
    myPlot = PlotlyJS.plot([PlotlyJS.scatter(x=mylist, y=valFH, name="valFH"),
                            PlotlyJS.scatter(x=mylist, y=valFVH, name="valFVH")])
    display(myPlot)

end

function influenceLW()
    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"
    pathWild = pathWD * "/WildLW/"
    myModelInput = pathWild * "/input.txt"

    paramName = "LW"
    paramList = [0.0, 10.0, 50.0, 100.0]
    readCreateNumericalModel(myModelInput, paramName, paramList)

    resultFolders = ["/home/hetier/Documents/PhD/SimulationNumerique/WildLW/wildRK4_", 
                     "/home/hetier/Documents/PhD/SimulationNumerique/WildLW/wildRK4_LW1",
                     "/home/hetier/Documents/PhD/SimulationNumerique/WildLW/wildRK4_LW2",
                     "/home/hetier/Documents/PhD/SimulationNumerique/WildLW/wildRK4_LW3",
                     "/home/hetier/Documents/PhD/SimulationNumerique/WildLW/wildRK4_LW4",
                     "/home/hetier/Documents/PhD/SimulationNumerique/WildLW/wildRK4_LW5"]

    resultFolders = solveWriteModels(pathWild)

    plotTrajectory1d([resultFolder * "/result.csv" for resultFolder in resultFolders],
                        pathWild*"/plot.html",
                        2;
                        title="",
                        legendTitle = latexstring("L_W \\text{ values:}"),
                        plotLegend = true,
                        legend=string.(append!([2.0], paramList)),
                        toPlot=true,
                        fontsize = 30)
end

function mainTest()
    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"
    pathAnthropo = pathWD * "/Anthropized/"

    myModel = readNumericalModel(pathAnthropo * "/input.txt")
    math = myModel.mathModel

    
end

# mainBifurcation()
mainSolveModel()
end
println(t, " seconds to execute the code")