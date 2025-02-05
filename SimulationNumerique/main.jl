t = @elapsed begin

abstract type mathematicalModel end
abstract type numericalModel end

include("modelHunterRC.jl")
include("modelHunterRCRK4.jl")
include("modelHunterRCNS.jl")
include("modelHunterRCNSImplicit.jl")
include("writer.jl")
include("reader.jl")

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
    pathHunter = pathWD * "/HunterRC/GradI"
    inputPath = pathHunter 
    outputPath = pathHunter
    path = pathHunter


    solveWriteModels(inputPath, outputPath)

    # solveWriteModel(inputPath * "/inputLC.txt", outputPath)

    # myModel = readNumericalModel(path * "/inputHF.txt")
    # println(computeLambdaMincI0(myModel.mathModel))
    # println(computeLambdaMaxcI0(myModel.mathModel))
    # println(longTermDynamic(myModel.mathModel))
    # println(thresholdFW(myModel.mathModel))
    # println(computeValMax(myModel.mathModel))
    # println(equilibriumHFW(myModel.mathModel))
    # println(computeDeltaStab(myModel.mathModel))
end

function mainBifurcation()
    pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique"
    pathHunter = pathWD * "/HunterRC/FromEEFtoEEHF/"
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

function mainLimitCylce()
    pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique"
    pathHunter = pathWD * "/HunterRC/stabcI0/DiagrammeEF"
    inputPath = pathHunter
    outputPath = pathHunter
    path = pathHunter

    # solveWriteModel(inputPath * "/input.txt", outputPath)

    myModel = readNumericalModel(path * "/input.txt")
    println(computeLambdaMincI0(myModel.mathModel))
    # println(computeLambdaMaxcI0(myModel.mathModel))
    println(longTermDynamic(myModel.mathModel))
    println(computeCharacteristicLC(myModel))
    # println(computeValMax(myModel.mathModel))
    # println(equilibriumHFW(myModel.mathModel))
    # println(computeDeltaStab(myModel.mathModel))
end


mainSolveModel()
# mainBifurcation()
# mainLambdaMinMax()
# mainLimitCylce()

end
println(t, " seconds to execute the code")