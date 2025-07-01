t = @elapsed begin
    allFileLocation = splitdir(splitdir(pwd())[1])[1]
   
    include(allFileLocation * "/modelHunterHT2RC.jl")
    include(allFileLocation * "/modelHunterHT2RCRK4.jl")
    include(allFileLocation * "/modelHunterHT2RCNS.jl")
    include(allFileLocation * "/writer.jl")
    include(allFileLocation * "/reader.jl")
    
    function solveWriteModel(inputFile::String, outputDir::String)
        """
        Solve the model corresponding to the inputFile
        Store the results (results.csv and copy of inputFile) on outputDir
        """
        myModel = readNumericalModel(inputFile)   
        
        solveModel(myModel)
        # println(longTermDynamic(myModel.mathModel))
        writeResult(myModel, outputDir)
    
        startInput = findlast("/input", inputFile)
        inputName = inputFile[first(startInput):end]
    
        # if outputDir * inputName != inputFile
        #     cp(inputFile, outputDir * inputName, force=true) 
        # end 
    end
    
    function solveWriteModels(inputDir::String, outputDirGeneral::String)
        """
        Search all the input files inputXXX.txt in InputDir, and solve the corresponding model
        Results (= result.csv and inputXXX.txt) are stored on outputDirGeneral/XXX/
        """
        searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
        inputNames = searchdir(inputDir, "input")
    
        for inputName in inputNames
            if isfile(inputDir * "/" * inputName)
                
                outputDir = outputDirGeneral * "/" * inputName[6:end-4]
                solveWriteModel(inputDir * "/" * inputName, outputDir)
            end
        end
    end
    
    function mainSolveModel()
        pathWD = pwd()
        inputPath = pathWD
        outputPath = pathWD
    
        inputFile = inputPath * "/inputLCI0.txt"
        model = readNumericalModel(inputFile)
        mathModel = model.mathModel
        println(equilibriumH(mathModel))
        println(equilibriumHFW(mathModel))
        
        println("N_{I=0} = ", thresholdN(mathModel))
        println("If lambda > lambda_{max, I=0}, convergence toward a limit cycle")
        println("lambda_{max, I=0} = ", computelambdaMaxI0(mathModel))
        println("lambda = ", mathModel.lambdaFWH)
        println(longTermDynamicStr(mathModel))

        solveWriteModel(inputFile, outputPath)
    
    end
    
mainSolveModel()


end
println("Code basé sur le fichier formulation4.tex")
println(t, " seconds to execute the code")