t = @elapsed begin

abstract type mathematicalModel end
abstract type numericalModel end

# include("modelWild.jl")
# include("modelWildRK4.jl")
# include("modelAnthropized.jl")
# include("modelAnthropizedRK4.jl")
# include("modelMigration.jl")
# include("modelMigrationRK4.jl")
include("modelHunterRC.jl")
include("modelHunterRCRK4.jl")
# include("modelTest.jl")
# include("modelTestRK4.jl")
# include("modelDomestic.jl")
# include("modelDomesticRK4.jl")

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
    pathHunter = pathWD * "/Hunter/"
    pathHunterTest = pathWD * "/HunterTest/"
    pathDomestic = pathWD * "/Domestic/"
    path = pathHunter

    myModel = readNumericalModel(path * "input.txt")
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

    dicEqNamesColor = Dict(["(F_W)" => "silver",
                            "(F_W)(H_DH_W)" => "firebrick",
                            "(F_W)(H_DF_WH_W)" => "green"])

    
    
    ## parameter list ##
    listLambdaFWH = [k for k = 0.001:0.05:0.5]
    listKF = [k for k=10.:1.:200]
    listLambdaVH = [k for k = 0.1:0.002:1]
    listLV = [k for k = 1.0:0.1:75.1]
    listbeta = [k for k = 0.1:0.01:12]
    lista = [k for k = 0.01:0.05:15]
    listb = [k for k = 0.01:0.05:15]
    listc = [k for k = 1.5:0.01:2.1]
    listDelta1 = [k for k = 0.01:0.001:0.5]

    ## dir and file name ##
    
    dirName = path * "/bifurcation/"
    nameSuffix = "bifurcationDiagramCLambdaFWHTest"
    println(math)

    ## computation ##
    bifurcationMatrix = computeBifurcationDiagram(math, "lambdaFWH", listLambdaFWH, "c", listc)
    filePrefix = writeDiagram(math, bifurcationMatrix, dirName; 
                                        nameSuffix = nameSuffix)
    ## plot ##
    data, layout = plotBifurcationFile(filePrefix * ".csv";
                        toSave = false, 
                        saveFileName = filePrefix * ".html",
                        # title = latexstring("\\text{Bifurcation Diagram in } (c, L_V) \\text{ plane}"),
                        title = "",
                        xlabel = latexstring("\\LARGE{\\lambda_{FH}}"),
                        ylabel = latexstring("\\LARGE{c}"),
                        dicEqColor = dicEqNamesColor,
                        toPlot = true,
                        eqVals = false,
                        titleEqVals = "\\text{: average of equilibrium values}",
                        fontsize = 30)

    # lambdamax = 4 * math.rF / (math.m * math.aW * math.KF)
    # println(lambdamax)
    # # addtraces(myPlot, scatter(;listc, lambdamax * ones(length(listc)), mode="lines"))
    # display(PlotlyJS.plot([data, PlotlyJS.scatter(x= listc, y = lambdamax * ones(length(listc)), mode="lines")], layout))              
end

function mainSolveModel()
    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"
    # pathWild = pathWD * "/Wild/"
    # pathAnthropo = pathWD * "/Anthropized/"
    # pathMigration = pathWD * "/Migration/"
    pathHunter = pathWD * "/HunterRC/"
    # pathTest = pathWD * "/Test/"
    # pathDomestic = pathWD * "/Domestic/"
    path = pathHunter
    
    # resultFolders = solveWriteModels(path)
    # println(resultFolders)
    resultFolder = solveWriteModel(path * "input.txt")
    resultFolders = [resultFolder]
    myModel = readNumericalModel(path * "input.txt")
    
    println("HDFWHW : ", equilibriumHDFW(myModel.mathModel))



    ## Plotting ##

    listColor = ["red", "red", "red", "red", "red", "red", "blue", "blue", "blue", "blue",
                "green", "green", "green"]
    # plotPhasePortrait([resultFolder * "/result.csv" for resultFolder in resultFolders],
    #                     path*"/plotOrbit.html",
    #                     [1,3];
    #                     # title=latexstring("\\text{Orbit in the } F_W - V_W \\text{ plane}"),
    #                     legend=Dict([1=>latexstring("EE^{F_AH_A}_1"), 7=>latexstring("EE^{F_AH_A}_1"), 11=>latexstring("TE")]),
    #                     listColor = listColor,
    #                     completeLegend = true,
    #                     showLegend = true,
    #                     toPlot=true,
    #                     fontsize = 30)
    trajectory, layout = plotPhasePortrait([resultFolder * "/result.csv" for resultFolder in resultFolders],
                        path * "/plotFVH.html",
                        [1,2];
                        # title=latexstring("\\text{Orbit in the } F_W - V_W \\text{ plane}"),
                        # legend=Dict([7=> "7", 1=>"1", 2=>"2",3=>"3",4=>"4",5=>"5", 6=>"6"]),
                        toPlot=true,
                        listColor = ["red",  "red", "red", "blue", "blue", "blue", "blue", "blue", "blue", "green"],
                        fontsize = 30, transform="plusa")
    # plotTrajectory1d([resultFolder * "/result.csv" for resultFolder in resultFolders],
    #                     path*"/plotHFV1D.html",
    #                     3;
    #                     # title=latexstring("\\text{Orbit in the } F_W - V_W \\text{ plane}"),
    #                     # legend=[latexstring("allee"), "ecoService"],
    #                     showLegend = true,
    #                     toPlot=true,
    #                     fontsize = 15)
end


function mainDeltaStab()
    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"
    pathHunter = pathWD * "/HunterRC/"
    path = pathHunter

    myModel = readNumericalModel(path * "input.txt")
    mModel = myModel.mathModel

    lambdaMax = (mModel.muD - mModel.rD) * mModel.rF / (mModel.m * mModel.c)
    println("LambdaMax: ", lambdaMax)

    println("DeltaStab : ", computeDeltaStab(mModel))

    listLambda = Vector(range(0.0001, lambdaMax, 1000))
    listAlpha = Vector(range(0.0,1.0,1000))

    # listDelta = [computeDeltaStab(mModel, alpha=alpha)[2] for alpha in listAlpha[1:end-1]]
    # mPlot = plot(listAlpha, listDelta)
    # display(mPlot)
    
    deltaStabMatrix, FeqMatrix = computeDeltaStabMatrix(mModel, listAlpha[1:end-1], listLambda[1:end-1])
    
    dirName = path * "/stab/"
    nameSuffixDelta = "stabDiagramLambdaRealTest"
    nameSuffixFeq = "FeqDiagramLambdaRealTest"

    filePrefix = writeDiagram(mModel, deltaStabMatrix, dirName; 
                nameSuffix = nameSuffixDelta)
    filePrefix = writeDiagram(mModel, FeqMatrix, dirName; 
                nameSuffix = nameSuffixFeq)

end

mainSolveModel()
# mainDeltaStab()
# mainBifurcation()

end
println(t, " seconds to execute the code")