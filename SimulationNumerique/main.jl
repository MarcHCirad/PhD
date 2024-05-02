t = @elapsed begin

include("modelFVHRK4.jl")
include("modelFVHNSS.jl")
include("modelEcoServiceRK4.jl")
include("modelEcoServiceFVRK4.jl")
include("modelAlleeEffect.jl")
include("modelAlleeEffectRK4.jl")
include("modelAlleeEffectPref.jl")
include("modelAlleeEffectPrefRK4.jl")
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
    # model = allModels["AlleeEffect"]
    # math = model.mathModel
    
    # table2 = modelTable(math; biologicalModel = true)
    # sort!(table2)
    # possibleeq = (unique(table2[!,:Equilibrium]))
    # println(size(possibleeq))
    # println(sort(possibleeq))
    # output = pathWDFV * "/"*string(typeof(model))*"/TableBio.csv"
    # CSV.write(output, table2)


    # listLambdaFH = [k for k = 0.001:0.01:0.5]
    # listLambdaVH = [k for k = 0.001:0.001:1]
    # listLV = [k for k = 20:0.1:70]
    # listbeta = [k for k = 0.01:0.01:7]

    dicEqNamesNbr = Dict(["(F)" =>"royalblue",
                        "(F)(FH)"=>"firebrick", 
                        "(F)(FH)(VH_2)"=>"darkorange",
                        "(F)(FH_β)"=>"green", "(F)(FH_β)(FH)"=>"orchid", "(F)(FH_β)(FH)(VH_2)"=>"black",
                        "(F)(FH_β)(VH_2)"=>"darkgreen", "(F)(FVH_β)"=>"silver", "(F)(FVH_β)(FH)"=>"chocolate",
                        "(F)(FVH_β)(FH)(VH_2)"=>"plum", "(F)(FVH_β)(VH_2)"=>"goldenrod", "(F)(H)"=>"gold",
                        "(F)(H)(VH_2)"=>"bisque", "(F)(VH_2)"=>"turquoise",
                        "(F)(VH_β)"=>"peru", "(F)(VH_β)(VH_2)"=>"chartreuse"])

    # bifurcationMatrix = computeBifurcationDiagram(math, "beta", listbeta, "LV", listLV)
    # dirName = writeBifurcationDiagram(model, bifurcationMatrix, pathWDFV)
    dirName = "/home/hetier/Documents/PhD/SimulationNumerique/Test/alleeEffectRK4"
    plotBifurcationFile(dirName * "/bifurcationDiagram.csv", 
                        pathWDFV*"/bifurcation.html" ;
                        dicEqColor = dicEqNamesNbr,
                        xlabel = "Beta",
                        ylabel = "L_V",
                        toPlot = true)
    

    # solveModel(allModels)
    # # F1, V1, H1 = eqFVH1(math)
    # # F2, V2, H2 = eqFVH2(math)
    # # println("F1 : ", F1, " V1 : ", V1, " H1 : ", H1)
    # # println("F2 : ", F2, " V2 : ", V2, " H2 : ", H2)
    # resultFolder = writeNumericalModel(allModels, pathWDFV)
    # println("Results have been writting in the following folders :")
    # for dir in resultFolder
    #     println(dir)
    # end
    
    # plotTrajectory1d([dir * "/result.csv" for dir in resultFolder],
    #                     pathWDFV*"/plot.html",
    #                     title="Comparison of two models : Eq FH",
    #                     legend=["allee", "ecoService"],
    #                     toPlot=true)

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