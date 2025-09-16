t = @elapsed begin

abstract type mathematicalModel end
abstract type numericalModel end

    allFileLocation = splitdir(splitdir(pwd())[1])[1]
   
    include(allFileLocation * "/modelHunterHT2RC.jl")
    include(allFileLocation * "/modelHunterHT2RCRK4.jl")
    include(allFileLocation * "/writer.jl")
    include(allFileLocation * "/reader.jl")

function mainBifurcation()
    pathWD = pwd() * "/"

    inputFile = pathWD * "inputDiagram.txt"
    myModel = readNumericalModel(inputFile)
    mathModel = myModel.mathModel
    println(equilibriumHFW(mathModel))

    ### 2D
    param1Name = "alpha"
    listAlpha = Vector(range(0.0,0.99,2000))
    param2Name = "lambdaFWH"
    listLambda = Vector(range(0.000001, 0.2, 5000))


    bifurcationDiagram = computeBifurcationDiagram(mathModel, param1Name,
                                                listAlpha[1:end], param2Name, listLambda[1:end])
    lambdaMin, lambdaMax = computeLambdaMatrixAlpha(mathModel, listAlpha)

    dirName = pathWD

    filePrefix = writeDiagram2D(dirName, inputFile, bifurcationDiagram ; nameSuffix = "Diagram")
    writeDiagram1D(dirName, inputFile, lambdaMin; nameSuffix = "Diagram1DMin")
    writeDiagram1D(dirName, inputFile, lambdaMax; nameSuffix = "Diagram1DMax")
    

end
  
    
mainBifurcation()
    
    
end
println("Code basé sur le fichier formulation4.tex")
println(t, " seconds to execute the code")