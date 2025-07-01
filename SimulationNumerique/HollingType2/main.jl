t = @elapsed begin

    abstract type mathematicalModel end
    abstract type numericalModel end
    
    include("/home/hetier/Documents/0PhD/SimulationNumerique/modelHunterHT2RC.jl")
    include("/home/hetier/Documents/0PhD/SimulationNumerique/modelHunterHT2RCRK4.jl")
    include("/home/hetier/Documents/0PhD/SimulationNumerique/writer.jl")
    include("/home/hetier/Documents/0PhD/SimulationNumerique/reader.jl")
    
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
    
        if outputDir * inputName != inputFile
            cp(inputFile, outputDir * inputName, force=true) 
        end 
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
        pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique/HollingType2/"
        inputPath = pathWD
        outputPath = pathWD * "/LCI0/"
    
        # solveWriteModels(inputPath, outputPath)
        inputFile = inputPath * "/inputLCI0.txt"
        model = readNumericalModel(inputFile)
        mathModel = model.mathModel
        println("lambda Max : ", computelambdaMaxI0(mathModel))
        solveWriteModel(inputFile, outputPath)
    
    end
    
mainSolveModel()


end
println("Code basé sur le fichier formulation4.tex")
println(t, " seconds to execute the code")