function writeResult(myModel::numericalModel, dirName::String)
    if !isdir(dirName)
        mkdir(dirName)
    end
    resultFile = dirName * "/result.csv"

    CSV.write(resultFile, 
        DataFrame(transpose(myModel.result), :auto),
        header = myModel.variablesNames)
end


function writeDiagram(dirPrefix::String,
                    inputFile::String,
                    diagram::Matrix{Any};
                    nameSuffix = "Diagram")
    
    if !isdir(dirPrefix)
        mkpath(dirPrefix)
    end


    indInput = first(findfirst("input", inputFile))
    copyInputFile = dirPrefix * "/" * inputFile[first(indInput):end]
    if copyInputFile != inputFile
        cp(inputFile, copyInputFile, force=true)
    end

    fileName = dirPrefix * "/" * nameSuffix * ".csv"
    println(fileName)

    header::Vector{String} = diagram[1, :]
    df = DataFrame(diagram[2:end, :], :auto)
    rename!(df, header)

    CSV.write(fileName, df)
    
    return fileName[1:end-4]
end
