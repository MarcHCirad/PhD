function writeMathematicalModel(myModel::MathematicalModel, dirName::String)
    if !isdir(dirName)
        mkdir(dirName)
    end

    parameterFile = dirName * "/math_parameters.txt"
    open(parameterFile, "w") do fout
        write(fout, string(typeof(myModel))*"\n")
        parameters = fieldnames(typeof(myModel))
        for ind in eachindex(parameters)
            parameter = parameters[ind]
            value = getfield(myModel, parameter)
            write(fout, string(parameter)*"\t"*string(value)*"\n")
        end
    end
end

function writeResult(myModel::numericalModel, dirName::String)
    if !isdir(dirName)
        mkdir(dirName)
    end
    resultFile = dirName * "/result.csv"

    CSV.write(resultFile, 
        DataFrame(transpose(myModel.result), :auto),
        header=["time", "F", "V", "H"])
end


function writeNumericalModel(myModel::numericalModel, dirName::String)
    if !isdir(dirName)
        mkdir(dirName)
    end
    numParamterFile = dirName * "/numerical_parameters.txt"
    parameters = fieldnames(typeof(myModel))

    mathModel = getfield(myModel, parameters[1])
    writeMathematicalModel(mathModel, dirName)

    open(numParamterFile, "w") do fout
        write(fout, string(typeof(myModel))*"\n")
        for ind in 2:length(parameters)-1
            parameter = parameters[ind]
            value = getfield(myModel, parameter)
            write(fout, string(parameter)*"\t"*string(value)*"\n")
        end
    
    end

    writeResult(myModel, dirName)
end

function createMathModel(modelType::String, param::Dict{String, Float64})
    if modelType == "FVH"
        nbrParameters = length(fieldnames(modelFVH))
        if nbrParameters != length(param)
            println("Number of parameter is incoherent with model")
            return
        end
        myModel = modelFVH(param)
        return myModel
    elseif modelType == "RV"
        nbrParameters = length(fieldnames(modelRV))
        if nbrParameters != length(param)
            println("Number of parameter is incoherent with model")
            return
        end
        myModel = modelRV(param)
        return myModel
    elseif modelType == "EcoService"
        nbrParameters = length(fieldnames(modelEcoService))
        if nbrParameters != length(param)
            println("Number of parameter is incoherent with model")
            return
        end
        myModel = modelEcoService(param)
        return myModel
    else
        println("Model type : ", modelType, " is not implemented.")
        return
    end
end

function createNumericalModel(numericalModelType::String, mathModelType::String, 
        mathModel::MathematicalModel, numericalParam::Dict{String, Float64}, 
        initialValues::Dict{String, Float64})

    if mathModelType == "FVH"
        if numericalModelType == "RK4"
            myModel = FVHRK4(mathModel, numericalParam, initialValues)
        elseif numericalModelType == "NSS"
            myModel = FVHNSS(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "RV"
        if numericalModelType == "RK4"
            myModel = RVRK4(mathModel, numericalParam, initialValues)
        elseif numericalModelType == "NSS"
            myModel = RVNSS(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "EcoService"
        if numericalModelType == "RK4"
            myModel = ecoServiceRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    end

end

function readNumericalModel(dirName::String)
    fileName = dirName * "test.txt"
    myCSV = CSV.read(fileName, DataFrame, header=false, delim=",", missingstring="NA")
    
    ind = 1
    myCSV[ind, 1] ## Should be ## Mathematical Parameters ##
    ind += 1
    mathModelType::String = myCSV[ind, 2]
    mathParam = Dict{String, Float64}()
    ind += 1
    while myCSV[ind, 1] != "## Initial Values ##"
        mathParam[myCSV[ind, 1]] = tryparse(Float64, myCSV[ind, 2])
        ind += 1
        println(myCSV[ind, 1])
    end
    mathModel = createMathModel(mathModelType, mathParam)

    ind += 1
    initialValues = Dict{String, Float64}()
    while myCSV[ind, 1] != "## Numerical Parameters ##"
        initialValues[myCSV[ind, 1]] = tryparse(Float64, myCSV[ind, 2])
        ind += 1
    end
    ind += 1

    myCSV[ind, 1] ## Should be Type
    numericalModelType::String = myCSV[ind, 2]
    numericalParam = Dict{String, Float64}()

    ind += 1

    while ind <= nrow(myCSV)
        numericalParam[myCSV[ind, 1]] = tryparse(Float64, myCSV[ind, 2])
        ind += 1
    end
    
    numericalModel = createNumericalModel(numericalModelType, mathModelType, 
                    mathModel, numericalParam, initialValues)
    
    return numericalModel
end


function readNumericalModelTest(dirName::String)
    fileName = dirName * "/input.txt"
    open(fileName, "r") do fin    

        readline(fin) ## Should be ## Mathematical Parameters ##

        line = readline(fin)
        indTab = last(findfirst("    ", line))
        mathModelType::String = line[last(indTab)+1:end]
        mathParam = Dict{String, Float64}()
        
        line = readline(fin)
        while line != "## Initial Values ##"
            indTab = findfirst("    ", line)
            mathParam[line[1:first(indTab)-1]] = tryparse(Float64, line[last(indTab)+1:end])
            line = readline(fin)
        end
        mathModel = createMathModel(mathModelType, mathParam)

        line = readline(fin)

        initialValues = Dict{String, Float64}()
        while line != "## Numerical Parameters ##"
            indTab = findfirst("    ", line)
            initialValues[line[1:first(indTab)-1]] = tryparse(Float64, line[last(indTab)+1:end])
            line = readline(fin)
        end
        
        line = readline(fin)
        indTab = last(findfirst("    ", line))
        numericalModelType::String = line[last(indTab)+1:end]
        numericalParam = Dict{String, Float64}()

        for line in eachline(fin)
            indTab = findfirst("    ", line)
            numericalParam[line[1:first(indTab)-1]] = tryparse(Float64, line[last(indTab)+1:end])
        end

        numericalModel = createNumericalModel(numericalModelType, mathModelType, 
                mathModel, numericalParam, initialValues)

        return numericalModel
    end
end
