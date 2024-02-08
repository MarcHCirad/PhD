t = @elapsed begin

include("modelFVHRK4.jl")
include("modelFVHNSS.jl")
include("modelRVRK4.jl")
include("modelRVNSS.jl")
include("modelEcoService.jl")
include("modelEcoServiceRK4.jl")
include("plotFile.jl")

function main()
    ## Model parameters
    modelParam = Dict([("rV",0.2), ("KV",50.), ("alpha",0.01), ("muV",0.1), 
        ("rF",0.91), ("KF",429.2), ("omega",0.1), ("f",0.2), ("muF",0.1)])
    modelParam["lambdaVH"] = 0.0332
    modelParam["lambdaFH"] = 0.001
    modelParam["gamma"] = 0.9
    modelParam["betaF"] = 5.
    modelParam["betaH"] = 1.
    modelParam["g"] = 0.03
    modelParam["a"] = 5.5
    modelParam["b"] = 1.
    modelParam["c"] = 10.
    modelParam["rH"] = 0.05
    
    ## Create a mathematical model
    myMathModel = modelEcoService(modelParam)
    println(interpretExistenceTresholds(myMathModel))

    # ## Numerical parameters
    t0, tf = 0., 500.
    n = 50000
    dt = (tf-t0)/n
    numericalParam = Dict([("t0", t0), ("tf", tf), ("dt", dt)])

    F0, V0, H0 = 50.2,0.1, 100.
    initialValues = Dict([("F0", F0), ("V0", V0), ("H0", H0)])

    # ## Create a numerical model 
    # myModel = ecoServiceRK4(modelParam, numericalParam, initialValues)
    # solveModel(myModel)
    # CSV.write("model1.csv", 
    #     DataFrame(transpose(myModel.result), :auto),
    #     header=["time", "F", "V", "H"])

    # modelParam["rH"] = 0.02
    # modelParam["g"] = 1
    # myModel2 = ecoServiceRK4(modelParam, numericalParam, initialValues)
    # solveModel(myModel2)
    # CSV.write("model2.csv", 
    #     DataFrame(transpose(myModel2.result), :auto),
    #     header=["time", "F", "V", "H"])

    
    
    plotTrajectory3d(["model1.csv", "model2.csv"], "test.html",
                    toPlot = true,
                    title="test",
                    legend=["rH = 0.05", "rH=0.02"])
end

main()

end
println(t, " seconds to execute the code")