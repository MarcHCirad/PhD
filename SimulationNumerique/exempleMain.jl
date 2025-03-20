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

function mainSolveModel()
    pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique"
    pathHunter = pathWD * "/PositiveFeedback/Orbit"
    inputPath = pathHunter 
    outputPath = pathHunter * "/LCI20"
    path = pathHunter


    # solveWriteModels(inputPath, outputPath)

    # solveWriteModel(inputPath * "/inputBeta0.txt", outputPath)

end

function mainBifurcation()
    pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique"
    pathHunter = pathWD * "/PositiveFeedback/"
    path = pathHunter

    inputFile = path * "inputTest.txt"
    myModel = readNumericalModel(inputFile)
    mModel = myModel.mathModel

    println(computeLambdaMax(mModel))
   
    listLambda = Vector(range(0.0001, 0.75, 2000))
    listAlpha = Vector(range(0.,0.98,1000))

    bifurcationDiagram = computeBifurcationDiagram(mModel, 
                                                listAlpha[1:end-1], listLambda[1:end-1])

    dirName = path * "BifurcationTest"

    filePrefix = writeDiagram(dirName, inputFile, bifurcationDiagram; nameSuffix="BifurcationTest")
    # filePrefix = writeDiagram(dirName, inputFile, lambdaMin; nameSuffix="LambdaMin")
    # filePrefix = writeDiagram(dirName, inputFile, lambdaMax; nameSuffix="LambdaMax")
end

function mainBifurcation1D()
    pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique"
    path = pathWD * "/Bifurcation1D"

    inputFile = path * "/inputBeta0.txt"
    myModel = readNumericalModel(inputFile)
    
    mathModel = myModel.mathModel
    paramName = "lambdaFWH"
    paramTable = Vector(range(0.00,0.0042,2000))

    bifurcationDiagram = bifurcationEquilibriumValues(mathModel, 
                                                paramName, paramTable)

    dirName = path * "/BifurcationTest"

    filePrefix = writeDiagram(dirName, inputFile, bifurcationDiagram)
end


function mainLambdaMinMax()
    pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique"
    pathHunter = pathWD * "/HunterRC/"
    path = pathHunter
    inputFile = path * "stabcI0/Plot3/input.txt"
    myModel = readNumericalModel(inputFile)
    mModel = myModel.mathModel

    listMW = Vector(range(0.10001,1.2,2000))
    listAlpha = Vector(range(0.,0.70,2000))
    
    LambdaExistence, LambdaStab = computeLambdaMatrixAlphaMW(mModel, 
                            listAlpha[1:end-1],
                            listMW[1:end-1])

    dirName = path * "stabcI0/Plot3"

    nameSuffixExistence = "lambdaMAlphaExistence"
    nameSuffixStab = "lambdaMAlphaStab"

    filePrefix = writeDiagram(dirName, inputFile, LambdaExistence; nameSuffix=nameSuffixExistence)
    filePrefix = writeDiagram(dirName, inputFile, LambdaStab; nameSuffix=nameSuffixStab)


end




# mainSolveModel()
mainBifurcation1D()
# mainBifurcation()
# mainLambdaMinMax()

end
println("Code basé sur le fichier formulation4.tex")
println(t, " seconds to execute the code")