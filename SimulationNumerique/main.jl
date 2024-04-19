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

    pathWD = "/home/hetier/Documents/PhD/SimulationNumerique"

    pathWDFV = pathWD * "/Test"
    allModels = readNumericalModel(pathWDFV)
    model = allModels["AlleeEffect"]
    math = model.mathModel
    
    # table2 = modelTable(math)
    # sort!(table2)
    # possibleeq = (unique(table2[!,:Equilibrium]))
    # println(size(possibleeq))
    # println(possibleeq)
    # output = pathWDFV * "/"*string(typeof(model))*"/TableFH.csv"
    # CSV.write(output, table2)

    dictThresholds = computeThresholds(math)
    println("Equilibrium : ", interpretThresholds(math, dictThresholds))
    # listLambdaFH = [k for k = 0.001:0.001:0.5]
    # listLambdaVH = [k for k = 0.001:0.001:0.5]
    # listc = [k for k = 0.001:0.01:2]
    # listbeta = [k for k = 0.001:0.01:2]
    # bifurcationMatrix = computeBifurcationDiagram(math, "c", listc, "beta", listbeta)
    # dirName = writeBifurcationDiagram(model, bifurcationMatrix, pathWDFV)
    # plotBifurcationFile(dirName * "/bifurcationDiagram.csv", 
    #                     pathWDFV*"/bifurcation.html" ;
    #                     xlabel = "c",
    #                     ylabel = "beta",
    #                     toPlot = true)
    

    solveModel(allModels)
    println(eqVH1(math))
    println(eqVH2(math))
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