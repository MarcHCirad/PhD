t = @elapsed begin
    allFileLocation = splitdir(splitdir(pwd())[1])[1]
   
    include(allFileLocation * "/modelHunterHT2RC.jl")
    include(allFileLocation * "/modelHunterHT2RCRK4.jl")
    include(allFileLocation * "/writer.jl")
    include(allFileLocation * "/reader.jl")
    
    function solveWriteModel(inputFile::String, outputDir::String; nameSuffix = "")
        """
        Solve the model corresponding to the inputFile
        Store the results (results.csv and copy of inputFile) on outputDir
        """

        myModel = readNumericalModel(inputFile)   
        
        solveModel(myModel)
        writeResult(myModel, outputDir; nameSuffix = nameSuffix)
    end
    
    function mainSolveModel()
        pathWD = pwd()
        inputPath = pathWD
        outputPath = pathWD

        suffix = "" #for inputLC.txt
        # suffix = "Int" #for inputLCInt.txt
        # suffix = "Ext" #for inputLCExt.txt

        inputFile = inputPath * "/inputLC" * suffix * ".txt"
        model = readNumericalModel(inputFile)
        mathModel = model.mathModel
        # println(equilibriumH(mathModel))
        println(equilibriumHFW(mathModel))
        println("Beta Max :", computeBetaMax(mathModel))
        
        if mathModel.I == 0
            println("N_{I=0} = ", thresholdNI0(mathModel))
            println("If lambda > lambda_{max, I=0}, convergence toward a limit cycle")
            println("lambda_{max, I=0} = ", computeLambdaMaxI0(mathModel))
            println("lambda = ", mathModel.lambdaFWH)
        end
        println(longTermDynamicStr(mathModel))

        solveWriteModel(inputFile, outputPath, nameSuffix = "LC" * suffix)
    
    end
    
mainSolveModel()


end
println("Code basé sur le fichier formulation4.tex")
println(t, " seconds to execute the code")