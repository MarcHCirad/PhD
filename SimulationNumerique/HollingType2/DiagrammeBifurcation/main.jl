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
            outputDir = outputDirGeneral * "/" * inputName[1:end-4]
            solveWriteModel(inputDir * "/" * inputName, outputDir)
        end
    end
end

function mainBifurcation()
    pathWD = pwd() * "/"

    listI = ["01"]#, "001", "01", "1"]

    for I in listI
        inputFile = pathWD * "inputDiagram" * string(I) * ".txt"
        myModel = readNumericalModel(inputFile)
        mathModel = myModel.mathModel
        # printAllCharacteristicValues(mathModel)
        println(equilibriumHFW(mathModel))

        ### 2D
        param1Name = "alpha"
        listAlpha = Vector(range(0.0,0.99,2000))
        param2Name = "lambdaFWH"
        listLambda = Vector(range(0.000001, 0.5, 5000))
        # listLambda = [0.0002261080540270135]
        # listAlpha = [0.01]
        # computeDeltaStab(mathModel)
        
        # eq = equilibriumHFW(mathModel)
        # println(eq)
        # computeDeltaStab3D(mathModel, eq[1])

        bifurcationDiagram = computeBifurcationDiagram(mathModel, param1Name,
                                                    listAlpha[1:end], param2Name, listLambda[1:end])
        lambdaMin, lambdaMax = computeLambdaMatrixAlpha(mathModel, listAlpha)

        dirName = pathWD

        filePrefix = writeDiagram2D(dirName, inputFile, bifurcationDiagram ; nameSuffix = "Diagram" * string(I))
        writeDiagram1D(dirName, inputFile, lambdaMin; nameSuffix = "Diagram1D" * string(I) * "Min")
        writeDiagram1D(dirName, inputFile, lambdaMax; nameSuffix = "Diagram1D" * string(I) * "Max")
    end

end

function mainSolveModel()
    pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique/PositiveFeedback/Bifurcation/"

    inputFile = pathWD * "input.txt"
    outputDir = pathWD * "Orbit/"

    myModel = readNumericalModel(inputFile)
    mathModel = myModel.mathModel
    printAllCharacteristicValues(mathModel)
    
    solveWriteModel(inputFile, outputDir)

end
    
    
mainBifurcation()
    
    
end
println("Code basé sur le fichier formulation4.tex")
println(t, " seconds to execute the code")