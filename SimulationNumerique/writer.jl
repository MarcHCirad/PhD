
function writeMathematicalModel(myModel::mathematicalModel, dirName::String)
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


function writeNumericalModel(myModel::numericalModel, dirPrefix::String)
    if !isdir(dirPrefix)
        mkdir(dirPrefix)
    end
    dirName = dirPrefix * "/"*string(typeof(myModel))
    numParamterFile =  dirName * "/numerical_parameters.txt"
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
    return dirName
end