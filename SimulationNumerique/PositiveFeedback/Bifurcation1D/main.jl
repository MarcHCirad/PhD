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

function mainBifurcation1D()
    pathWD = "/home/hetier/Documents/0PhD/SimulationNumerique/Bifurcation1D"

    inputFile = pathWD * "/inputBeta0.txt"
    myModel = readNumericalModel(inputFile)
    mathModel = myModel.mathModel
    
    paramName = "lambdaFWH"
    paramTable = Vector(range(0.00,0.0042,2000))

    bifurcationDiagram = bifurcationEquilibriumValues(mathModel, 
                                                paramName, paramTable)

    dirName = path * "/BifurcationTest"

    filePrefix = writeDiagram(dirName, inputFile, bifurcationDiagram)
end

mainBifurcation1D()


end
println("Code basé sur le fichier formulation4.tex")
println(t, " seconds to execute the code")