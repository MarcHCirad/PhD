t = @elapsed begin

include("modelFVHRK4.jl")
include("modelFVHNSS.jl")
include("modelRVRK4.jl")
include("modelRVNSS.jl")
include("modelEcoService.jl")
include("plotFile.jl")

function main()
    ## Model parameters
    modelParam = Dict([("rV",1.), ("KV",50.), ("alpha",0.02), ("muV",0.1), ("rF",0.71), ("KF",429.2), ("omega",0.25), ("f",0.25), ("muF",0.1), ("e",0.5)])
    modelParam["muH"] = 0.2
    modelParam["lambdaVH"] = 0.0332
    modelParam["lambdaFH"] = 0.01
    modelParam["H0"] = 10
    modelParam["gamma"] = 0.5
    modelParam["betaF"] = 1.
    modelParam["betaH"] = 1.
    modelParam["g"] = 1.
    modelParam["a"] = 0.5
    modelParam["b"] = 2.
    modelParam["c"] = 3.
    modelParam["rH"] = 2.5
    
    ## Create a mathematical model
    myMathModel = modelEcoService(modelParam)
    println(interpretExistenceTresholds(myMathModel))

    # start = 0.001
    # step = 0.001
    # stop = 0.2
    # listLambdaVH = Vector(range(start, stop, step=step))
    # start = 0.001
    # step = 0.001
    # stop = 0.2
    # listLambdaFH = Vector(range(start, stop, step=step))

    ## Compute a bifurcation diagram and write it in a csv file
    # computeBifurcationDiagram(myMathModel, listLambdaFH, listLambdaVH, "bifurcationDiagram.csv")
    ## Plot a bifurcation diagram from a csv file
    # plotBifurcationFile("bifurcationDiagram.csv", "bifurcationDiagram.html"; 
    #         xlabel = L"$\lambda_{FH}$",
    #         # ylabel = L"$\lambda_{VH}$",
    #         title = "Bifurcation Diagram",    
    #         toPlot=true)

    # ## Numerical parameters
    # t0, tf = 0., 500.
    # n = 50000
    # dt = (tf-t0)/n
    # numericalParam = Dict([("t0", t0), ("tf", tf), ("dt", dt)])

    # F0, V0, H0 = 8.,10.,5.
    # initialValues = Dict([("F0", F0), ("V0", V0), ("H0", H0)])

    # ## Create a numerical model 
    # myModel = RVNSS(modelParam, numericalParam, initialValues)
    # solveModel(myModel)
    # CSV.write("test.csv", 
    #     DataFrame(transpose(myModel.result), :auto),
    #     header=["time", "F", "V", "H"])
    # plotTrajectory("test.csv", "test.html",
    #                 toPlot = true,
    #                 title="test",
    #                 plotType="3d")
end

main()

end
println(t, " seconds to execute the code")