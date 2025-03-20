function createMathModel(modelType::String, param::Dict{String, Any})
    ## Check if the math model is implemented, and if it can be created. If yes, return it.
    if modelType == "Hunter"
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
    if mathModelType == "Hunter"
        if numericalModelType == "RK4"
            myModel = hunterRK4(mathModel, numericalParam, initialValues)
        elseif numericalModelType == "NS"
            myModel = hunterNS(mathModel, numericalParam, initialValues)
        elseif numericalModelType == "NSImplicit"
            myModel = hunterNSImplicit(mathModel, numericalParam, initialValues)
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
        mathParam = Dict{String, Any}()
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
