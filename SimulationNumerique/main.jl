t = @elapsed begin

include("modelFVHRK4.jl")
include("modelFVHNSS.jl")
include("modelRVRK4.jl")
include("modelRVNSS.jl")
include("modelEcoServiceRK4.jl")
include("modelEcoServiceFVRK4.jl")
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

    pathWDFV = pathWD * "/ModelEcoService/EquilibreEndemique"
    allModels = readNumericalModel(pathWDFV)
    solveModel(allModels)

    resultFolder = writeNumericalModel(allModels, pathWDFV)
    println("Results have been writting in the following folders :")
    for dir in resultFolder
        println(dir)
    end
    
    plotTrajectory1d([dir * "/result.csv" for dir in resultFolder],
                        pathWDFV*"/plot.html",
                        title="Comparison of two models : Eq FVH",
                        legend=["F,V, aF+bV,c", "F,V,H"],
                        toPlot=true)
end

main()

end
println(t, " seconds to execute the code")