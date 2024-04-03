t = @elapsed begin

include("modelFVHRK4.jl")
include("modelFVHNSS.jl")
include("modelEcoServiceRK4.jl")
include("modelEcoServiceFVRK4.jl")
include("modelAlleeEffect.jl")
include("modelAlleeEffectRK4.jl")
include("plotFile.jl")
include("writer.jl")
include("reader.jl")

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

function main()

    pathWD = "/home/hetier/Documents/Code/PhDCode/SimulationNumerique"

    pathWDFV = pathWD * "/Test"
    allModels = readNumericalModel(pathWDFV)    
    solveModel(allModels)

    resultFolder = writeNumericalModel(allModels, pathWDFV)
    println("Results have been writting in the following folders :")
    for dir in resultFolder
        println(dir)
    end
    
    plotTrajectory1d([dir * "/result.csv" for dir in resultFolder],
                        pathWDFV*"/plot.html",
                        title="Comparison of two models : Eq FH",
                        legend=["allee", "ecoService"],
                        toPlot=true)

    # plotPhasePortrait([dir * "/result.csv" for dir in resultFolder],
    #                     pathWDFV*"/plot.html",
    #                     [2,3];
    #                     title="Comparison of two models : Eq FH",
    #                     legend=["allee", "ecoService"],
    #                     toPlot=true)
end

main()

end
println(t, " seconds to execute the code")