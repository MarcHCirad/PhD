t = @elapsed begin

abstract type mathematicalModel end
abstract type numericalModel end

include("/home/hetier/Documents/0PhD/SimulationNumerique/modelHunterRC.jl")
include("/home/hetier/Documents/0PhD/SimulationNumerique/modelHunterRCRK4.jl")
include("/home/hetier/Documents/0PhD/SimulationNumerique/modelHunterRCNSImplicit.jl")
include("/home/hetier/Documents/0PhD/SimulationNumerique/writer.jl")
include("/home/hetier/Documents/0PhD/SimulationNumerique/reader.jl")

function solveWriteModel(inputFile::String, outputDir::String)
    myModel = readNumericalModel(inputFile)
    
    
    solveModel(myModel)
    println(longTermDynamic(myModel.mathModel))
    writeResult(myModel, outputDir)

    startInput = findlast("/input", inputFile)
    inputName = inputFile[first(startInput):end]

    if outputDir * inputName != inputFile
        cp(inputFile, outputDir * inputName, force=true) 
    end 
end

function solveWriteModels(inputDir::String, outputDirGeneral::String)
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
    pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique/PositiveFeedback/Bifurcation/"

    inputFile = pathWD * "input.txt"
    myModel = readNumericalModel(inputFile)
    mathModel = myModel.mathModel

    printAllCharacteristicValues(mathModel)

    ### 2D
    listLambda = Vector(range(0.0001, 2, 2000))
    listAlpha = Vector(range(0.,0.99,1000))

    bifurcationDiagram = computeBifurcationDiagram(mathModel, 
                                                listAlpha[1:end-1], listLambda[1:end-1])

    dirName = pathWD * "Diagram2D/"

    filePrefix = writeDiagram2D(dirName, inputFile, bifurcationDiagram)

    ### 1D
    paramName = "lambdaFWH"
    paramTable = Vector(range(0.0001,7,5000))

    bifurcationDiagram = bifurcationEquilibriumValues(mathModel, 
                                                paramName, paramTable)

    dirName = pathWD * "Diagram1D/"

    filePrefix = writeDiagram(dirName, inputFile, bifurcationDiagram)
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