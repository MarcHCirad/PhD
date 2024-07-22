t = @elapsed begin

abstract type mathematicalModel end
abstract type numericalModel end

include("modelWild.jl")
include("modelWildRK4.jl")
include("modelAnthropized.jl")
include("modelAnthropizedRK4.jl")

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

function mainBifurcation()

    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"

    pathWD = pathWD * "/Test2"
    allModels = readNumericalModel(pathWD)
    model = allModels["AlleeEffect"]
    math = model.mathModel   
    
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
                

    
    
    ## parameter list ##
    listLambdaFH = [k for k = 0.001:0.01:0.5]
    listLambdaVH = [k for k = 0.001:0.001:1]
    listLV = [k for k = 0.1:0.1:70]
    listbeta = [k for k = 0.1:0.01:12]
    lista = [k for k = 0.01:0.05:15]
    listb = [k for k = 0.01:0.05:15]
    listc = [k for k = 0.01:0.05:15]

    ## dir and file name ##
    dirName = pathWD * "/bifurcationAlleeEffect/CLV_test/"
    nameSuffix = "bifurcationDiagramCLV"
    filePrefix = dirName * "/" * nameSuffix

    ## computation ##
    bifurcationMatrix, bifF, bifV, bifH = computeBifurcationDiagram(math, "c", listc, "LV", listLV ; eqValues=true)
    filePrefix = writeBifurcationDiagram(math, bifurcationMatrix, dirName; 
                                        nameSuffix = nameSuffix,
                                        eqVals = true, bifF = bifF, bifV = bifV, bifH = bifH)

    ## plot ##
    plotBifurcationFile(filePrefix * ".csv", 
                        filePrefix * ".html" ;
                        # title = latexstring("\\text{Bifurcation Diagram in } (c, L_V) \\text{ plane}"),
                        title = "",
                        xlabel = latexstring("\\LARGE{c}"),
                        ylabel = latexstring("\\LARGE{L_V}"),
                        dicEqColor = dicEqNamesColor,
                        toPlot = true,
                        eqVals = false,
                        titleEqVals = "\\text{: average of equilibrium values}",
                        fontsize = 30)
end



function mainSolveModel()
    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"
    pathWild = pathWD * "/Wild/"
    
    resultFolders = ["/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_",
                    "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_2", 
                    "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_3", 
                    "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_4", 
                    "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_5",
                    "/home/hetier/Documents/PhD/SimulationNumerique/Wild/wildRK4_6"]
    resultFolders = solveWriteModels(pathWild)
    # myModel = readNumericalModel(pathWild * "/input.txt")
    # modelWild = myModel.mathModel
    # println(equilibriumFV(modelWild))

    ## Plotting ##
    # plotTrajectory1d([dir * "/result.csv" for dir in resultFolder],
    #                     nameFile,
    #                     title="Perturbation from equilibrium FVH : L_V decreased",
    #                     plotLegend = false,
    #                     legend=[latexstring("a_\\alpha")],
    #                     toPlot=true)

    plotPhasePortrait([resultFolder * "/result.csv" for resultFolder in resultFolders],
                        pathWild*"/plot.html",
                        [1,2];
                        title=latexstring("\\text{Orbit in the } F_W - V_W \\text{ plane}"),
                        # legend=[latexstring("allee"), "ecoService"],
                        toPlot=true)
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

mainSolveModel()

end
println(t, " seconds to execute the code")