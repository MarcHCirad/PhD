function createMathModel(modelType::String, param::Dict{String, Float64})
    ## Check if the math model is implemented, and if it can be created. If yes, return it.
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
        if nbrParameters > length(param)
            println("Number of parameter is incoherent with model EcoService")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelEcoService))
            return
        end
        myModel = modelEcoService(param)
        return myModel
    elseif modelType == "EcoServiceFV"
        nbrParameters = length(fieldnames(modelEcoServiceFV))
        if nbrParameters > length(param)
            println("Number of parameter is incoherent with model EcoServiceFV")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelEcoService))
            return
        end
        myModel = modelEcoServiceFV(param)
        return myModel
    elseif modelType == "AlleeEffect"
        nbrParameters = length(fieldnames(modelAlleeEffect))
        if nbrParameters - 1 > length(param)
            println("Number of parameter is incoherent with model AlleeEffect")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelAlleeEffect))
            return
        end
        myModel = modelAlleeEffect(param)
        return myModel
    elseif modelType == "AlleeEffectPref"
        nbrParameters = length(fieldnames(modelAlleeEffectPref))
        if nbrParameters - 1 > length(param)
            println("Number of parameter is incoherent with model AlleeEffectPref")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelAlleeEffectPref))
            return
        end
        myModel = modelAlleeEffectPref(param)
        return myModel
    elseif modelType == "WildV1"
        nbrParameters = length(fieldnames(modelWildV1))
        if nbrParameters > length(param)
            println("Number of parameter is incoherent with model WildV1")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelWildV1))
            return
        end
        myModel = modelWildV1(param)
        return myModel
    elseif modelType == "Wild"
        nbrParameters = length(fieldnames(modelWild))
        if nbrParameters -1 > length(param)
            println("Number of parameter is incoherent with model Wild")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelWild))
            return
        end
        myModel = modelWild(param)
        return myModel
    elseif modelType == "Anthropized"
        nbrParameters = length(fieldnames(modelAnthropized))
        if nbrParameters - 1 > length(param)
            println("Number of parameter is incoherent with model Anthropized")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelAnthropized))
            return
        end
        myModel = modelAnthropized(param)
        return myModel
    elseif modelType == "Migration"
        nbrParameters = length(fieldnames(modelMigration))
        if nbrParameters - 1 > length(param)
            println("Number of parameter is incoherent with model Migration")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelMigration))
            return
        end
        myModel = modelMigration(param)
        return myModel
    elseif modelType == "Hunter"
        nbrParameters = length(fieldnames(modelHunter))
        if nbrParameters - 3 > length(param)
            println("Number of parameter is incoherent with model Hunter")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelHunter))
            return
        end
        myModel = modelHunter(param)
        return myModel
    elseif modelType == "Domestic"
        nbrParameters = length(fieldnames(modelDomestic))
        if nbrParameters - 2 > length(param)
            println("Number of parameter is incoherent with model Domestic")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelDomestic))
            return
        end
        myModel = modelDomestic(param)
        return myModel
    elseif modelType == "Test"
        nbrParameters = length(fieldnames(modelTest))
        if nbrParameters - 1 > length(param)
            println("Number of parameter is incoherent with model Test")
            println("Numbers of parameters as arguments : ", length(param), " while expected number : ",
                        nbrParameters)
            println(keys(param), fieldnames(modelTest))
            return
        end
        myModel = modelTest(param)
        return myModel
    else
        println("Model type : ", modelType, " is not implemented.")
        return
    end
end

function createNumericalModel(numericalModelType::String, mathModelType::String, 
        mathModel::mathematicalModel, numericalParam::Dict{String, Float64}, 
        initialValues::Dict{String, Float64})
    ## Check if numerical model have been implemented for math model
    ## If yes, create the numerical model
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
    elseif mathModelType == "EcoServiceFV"
        if numericalModelType == "RK4"
            myModel = ecoServiceFVRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "AlleeEffect"
        if numericalModelType == "RK4"
            myModel = alleeEffectRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "AlleeEffectPref"
        if numericalModelType == "RK4"
            myModel = alleeEffectPrefRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "WildV1"
        if numericalModelType == "RK4"
            myModel = wildV1RK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "Wild"
        if numericalModelType == "RK4"
            myModel = wildRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "Anthropized"
        if numericalModelType == "RK4"
            myModel = anthropizedRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "Migration"
        if numericalModelType == "RK4"
            myModel = migrationRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "Hunter"
        if numericalModelType == "RK4"
            myModel = hunterRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "Domestic"
        if numericalModelType == "RK4"
            myModel = domesticRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    elseif mathModelType == "Test"
        if numericalModelType == "RK4"
            myModel = testRK4(mathModel, numericalParam, initialValues)
        else
            println("This type of numerical scheme is not implemented for this math. model")
        end
    end
end

function readNumericalModel(inputFile::String)
    """
    Take an directory name as input, conaining an input.txt file with appropriate format
    Return a list of numerical models, describe by the input.txt file
    """
    open(inputFile, "r") do fin
        readline(fin) ## Should be ## Mathematical Parameters ##

        ## Read the math model type(s)
        line = readline(fin)
        indTab = last(findfirst("    ", line)) ## escape the spaces
        mathModelType::String = line[last(indTab)+1:end]
      
        ## Read the math model parameters
        mathParam = Dict{String, Float64}()
        line = readline(fin)
        while line != "## Initial Values ##"
            indTab = findfirst("    ", line)
            mathParam[line[1:first(indTab)-1]] = tryparse(Float64, line[last(indTab)+1:end])
            line = readline(fin)
        end

        mathModel = createMathModel(mathModelType, mathParam)
        
        ## Read the initial values
        line = readline(fin)
        initialValues = Dict{String, Float64}()
        while line != "## Numerical Parameters ##"
            indTab = findfirst("    ", line)
            initialValues[line[1:first(indTab)-1]] = tryparse(Float64, line[last(indTab)+1:end])
            line = readline(fin)
        end
        
        ## Find the numerical model type
        line = readline(fin)
        indTab = last(findfirst("    ", line))
        numericalModelType::String = line[last(indTab)+1:end]
        numericalParam = Dict{String, Float64}()

        ## Read the numerical model parameters
        for line in eachline(fin)
            indTab = findfirst("    ", line)
            numericalParam[line[1:first(indTab)-1]] = tryparse(Float64, line[last(indTab)+1:end])
        end


        myNumModel = createNumericalModel(numericalModelType, mathModelType, 
                            mathModel, numericalParam, initialValues)
        return myNumModel
    end
end


function readMathModel(dirName::String)
    fileName = dirName * "/input.txt"
    open(fileName, "r") do fin    

        readline(fin) ## Should be ## Mathematical Parameters ##

        ## Read the math model type(s)
        line = readline(fin)
        indTab = last(findfirst("    ", line))
        mathModelTypesString::String = line[last(indTab)+1:end]

        commaIndices = findall(",", mathModelTypesString)
        firstIndices = append!([1], [last(ind)+1 for ind in commaIndices])
        lastIndices = append!([first(ind)-1 for ind in commaIndices], [length(mathModelTypesString)])
        mathModelTypes = [mathModelTypesString[firstIndices[k]:lastIndices[k]] for k in 1:(length(firstIndices))]
        
        ## Read the math model parameters
        mathParam = Dict{String, Float64}()
        line = readline(fin)
        while line != "## Initial Values ##"
            indTab = findfirst("    ", line)
            mathParam[line[1:first(indTab)-1]] = tryparse(Float64, line[last(indTab)+1:end])
            line = readline(fin)
        end

        mathModels = Dict{String, T where T<:mathematicalModel}()
        for mathModelType in mathModelTypes
            myModel = createMathModel(mathModelType, mathParam)
            mathModels[mathModelType] = myModel
        end
        return mathModels
    end
end