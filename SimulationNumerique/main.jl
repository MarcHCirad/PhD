t = @elapsed begin

include("plotFile.jl")
include("writer.jl")
include("reader.jl")

include("modelWildV1.jl")
include("modelWildV1RK4.jl")
include("modelWildV2.jl")
include("modelWildV2RK4.jl")

function solveModel(allModels::Dict{String, numericalModel})
    for model in collect(values(allModels)) 
        solveModel(model)
    end
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
    pathWDFV = pathWD * "/Wild"
    allModels = readNumericalModel(pathWDFV)
    resultFolder = ["/home/hetier/Documents/PhD/SimulationNumerique/Wild/RK4"]
    nameFile = pathWDFV*"/plotLVd.html"

    ## Solving ##
    solveModel(allModels)
    resultFolder = writeNumericalModel(allModels, pathWDFV)
    println("Results have been writting in the following folders :")
    for dir in resultFolder
        println(dir)
    end

    ## Plotting ##
    # plotTrajectory1d([dir * "/result.csv" for dir in resultFolder],
    #                     nameFile,
    #                     title="Perturbation from equilibrium FVH : L_V decreased",
    #                     plotLegend = false,
    #                     legend=[latexstring("a_\\alpha")],
    #                     toPlot=true)

    plotPhasePortrait([dir * "/result.csv" for dir in resultFolder],
                        pathWDFV*"/plot.html",
                        [1,2];
                        title="Comparison of two models : Eq FH",
                        legend=[latexstring("allee"), "ecoService"],
                        toPlot=true)
end

mainSolveModel()

end
println(t, " seconds to execute the code")