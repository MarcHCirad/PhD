function plotTrajectory1d(inputFileNames::Vector{String}, saveFileName::String;
        title::String = "",
        legend=Vector{String}(),
        toPlot = false)

    listColor = ["royalblue", "firebrick", "darkorange", "green", "orchid", "black"]

    if length(legend) < length(inputFileNames)
        append!(legend, ["File " * string(ind) for ind in length(legend):length(inputFileNames)])
    end

    fileName = inputFileNames[1]
    myCSV = CSV.read(fileName, DataFrame)
    header = names(myCSV)

    nameTime = header[1]
    nameColumn2 = header[2]
    nameColumn3 = header[3]
    nameColumn4 = header[4]

    time = myCSV[1:end, 1]
    column2 = myCSV[2:end, 2]
    column3 = myCSV[2:end, 3]
    column4 = myCSV[2:end, 4]

    lineCSV = attr(color=listColor[1])
    s1 = [PlotlyJS.scatter(x=time, y = column2, name=legend[1], line=lineCSV)]
    s2 = [PlotlyJS.scatter(x=time, y = column3, showlegend=false, line=lineCSV)]
    s3 = [PlotlyJS.scatter(x=time, y = column4, showlegend=false, line=lineCSV)]

    for ind in 2:length(inputFileNames)
        fileName = inputFileNames[ind]
        myCSV = CSV.read(fileName, DataFrame)
        lineCSV = attr(color=listColor[ind])
        time = myCSV[1:end, 1]

        column2 = myCSV[2:end, 2]
        push!(s1, PlotlyJS.scatter(x=time, y = column2, name=legend[ind], line=lineCSV))

        column3 = myCSV[2:end, 3]
        push!(s2, PlotlyJS.scatter(x=time, y = column3, showlegend=false, line=lineCSV))

        column4 = myCSV[2:end, 4]
        push!(s3, PlotlyJS.scatter(x=time, y = column4, showlegend=false, line=lineCSV))
    end


    l1 = PlotlyJS.Layout(yaxis_title = nameColumn2)
    l2 = PlotlyJS.Layout(yaxis_title = nameColumn3)
    l3 = PlotlyJS.Layout(xaxis_title = nameTime, yaxis_title = nameColumn4)
    p1 = PlotlyJS.plot(s1, l1)
    p2 = PlotlyJS.plot(s2, l2)
    p3 = PlotlyJS.plot(s3, l3)

    myPlot = [p1 ; p2 ; p3]
    PlotlyJS.relayout!(myPlot, title_text=title)
    if toPlot
        display(myPlot)
    end
        
    PlotlyJS.savefig(myPlot, saveFileName)
end

function plotTrajectory3d(inputFileNames::Vector{String}, saveFileName::String;
    variables = [],
    title = "",
    legend=Vector{String}(),
    toPlot = false)

    listColor = ["royalblue", "firebrick", "darkorange", "green", "orchid", "black"]

    if length(legend) < length(inputFileNames)
        append!(legend, ["File " * string(ind) for ind in length(legend):length(inputFileNames)])
    end

    fileName = inputFileNames[1]
    myCSV = CSV.read(fileName, DataFrame)
    header = names(myCSV)

    nameColumn2 = header[2]
    nameColumn3 = header[3]
    nameColumn4 = header[4]

    column2 = myCSV[2:end, 2]
    column3 = myCSV[2:end, 3]
    column4 = myCSV[2:end, 4]

    lineCSV = attr(color=listColor[1], width=2)
    markerCSV = attr(color=listColor[1], size=5, symbol="cross")
    layout = PlotlyJS.Layout(scene = attr(xaxis_title = nameColumn2, yaxis_title=nameColumn3, zaxis_title=nameColumn4), title=title)
    trajectory = [PlotlyJS.scatter(x=column2, y=column3, z=column4, type="scatter3d", name=legend[1], mode="lines", line=lineCSV)]
    push!(trajectory, PlotlyJS.scatter(x=[column2[1]], y=[column3[1]], z=[column4[1]], 
                type="scatter3d", mode="markers", marker=markerCSV, showlegend=false))
    

    for ind in 2:length(inputFileNames)
        fileName = inputFileNames[ind]
        myCSV = CSV.read(fileName, DataFrame)
        lineCSV = attr(color=listColor[ind], width=2)
        markerCSV = attr(color=listColor[ind], size=5, symbol="cross")

        column2 = myCSV[2:end, 2]
        column3 = myCSV[2:end, 3]
        column4 = myCSV[2:end, 4]
        push!(trajectory, PlotlyJS.scatter(x=column2, y=column3, z=column4, type="scatter3d", name=legend[ind], mode="lines", line=lineCSV))
        push!(trajectory, PlotlyJS.scatter(x=[column2[1]], y=[column3[1]], z=[column4[1]], 
                                                    type="scatter3d", mode="markers", marker=markerCSV, showlegend=false))
    end

    myPlot = PlotlyJS.plot(trajectory, layout)
    if toPlot
        display(myPlot)
    end


    PlotlyJS.savefig(myPlot, saveFileName)
end

function plotBifurcationFile(inputFileName::String, saveFileName::String;
        dicEqColor = Dict{String, String}(),
        xlabel = "",
        ylabel = "",
        title = "",    
        toPlot=false)
    """
    Take a CSV bifurcation file as input and produce the corresponding plot
    """
    myCSV = CSV.read(inputFileName, DataFrame)
    
    listParamX = tryparse.(Float64, myCSV[2:end, 1])
    listParamY = tryparse.(Float64, Vector(myCSV[1, 2:end]))

    eqNames = Matrix(myCSV[2:end, 2:end])

    eqNamesUnique = sort(unique(eqNames)) # Names of the different stability
    legendSize = length(eqNamesUnique)  # Number of the different stability

    dicEqNbr = Dict([eqNamesUnique[i+1]=>i for i in 0:legendSize-1])
    dicNbrEq = Dict([i => eqNamesUnique[i+1] for i in 0:legendSize-1])
    color = transpose(similar(eqNames, Float64))   # Matrix of color

    for ind_col in 1:length(listParamY)
        for ind_row in 1:length(listParamX)
            color[ind_col, ind_row] = dicEqNbr[eqNames[ind_row, ind_col]] # Fill the matrix
        end
    end

    ## Plot configuration
    colorStep = 1/legendSize
    myTickvals = [(2*k-1) / 2 for k = 1:legendSize]
    myTicktext = sort(collect(keys(dicEqNbr)))
    listColor = ["seashell", "darkred", "seagreen"]
    
    myColorScale = []
    counter = 1
    for k in 1:legendSize
        if dicNbrEq[k-1] in keys(dicEqColor)
            push!(myColorScale, ((k-1)*colorStep, dicEqColor[dicNbrEq[k-1]]))
            push!(myColorScale, ((k)*colorStep, dicEqColor[dicNbrEq[k-1]]))
        else
            push!(myColorScale, ((k-1)*colorStep, listColor[counter]))
            push!(myColorScale, ((k)*colorStep, listColor[counter]))
            counter += 1
        end

    end
    println(values(dicEqNbr))
    myzmin = -0.001
    myzmax = legendSize
    
    
    layout = PlotlyJS.Layout(xaxis_title= xlabel, yaxis_title=ylabel, title=title)
    data = PlotlyJS.heatmap(x = listParamX, y = listParamY, z = color,
                            colorbar=attr(tickmode="array",
                                            tickvals=myTickvals,
                                            ticktext=myTicktext),
                            autocolorscale = false,
                            colorscale = myColorScale,
                            zauto = false,
                            zmin = myzmin,
                            zmax = myzmax)

    myPlot = PlotlyJS.plot(data, layout)

    PlotlyJS.savefig(myPlot, saveFileName)

    if toPlot
        display(myPlot)
    end

end

function plotPhasePortrait(inputFileNames::Vector{String}, saveFileName::String,
    variables::Vector{Int64};
    title = "",
    legend=Vector{String}(),
    toPlot = false)

    listColor = ["royalblue", "firebrick", "darkorange", "green", "orchid", "black"]

    if length(legend) < length(inputFileNames)
        append!(legend, ["File " * string(ind) for ind in length(legend):length(inputFileNames)])
    end

    fileName = inputFileNames[1]
    myCSV = CSV.read(fileName, DataFrame)
    header = names(myCSV)

    indColumn2 = variables[1]+1
    indColumn3 = variables[2]+1

    nameColumn2 = header[indColumn2]
    nameColumn3 = header[indColumn3]

    column2 = myCSV[2:end, indColumn2]
    column3 = myCSV[2:end, indColumn3]

    lineCSV = attr(color=listColor[1], width=2)
    markerCSV = attr(color=listColor[1], size=5, symbol="cross")
    
    
    
    trajectory = [PlotlyJS.scatter(x=column2, y=column3, type="scatter1d", name=legend[1], mode="lines", line=lineCSV)]
    push!(trajectory, PlotlyJS.scatter(x=[column2[1]], y=[column3[1]], 
                type="scatter1d", mode="markers", marker=markerCSV, showlegend=false))
    

    for ind in 2:length(inputFileNames)
        fileName = inputFileNames[ind]
        myCSV = CSV.read(fileName, DataFrame)
        lineCSV = attr(color=listColor[ind], width=2)
        markerCSV = attr(color=listColor[ind], size=5, symbol="cross")

        column2 = myCSV[2:end, indColumn2]
        column3 = myCSV[2:end, indColumn3]
        push!(trajectory, PlotlyJS.scatter(x=column2, y=column3, type="scatter1d", name=legend[ind], mode="lines", line=lineCSV))
        push!(trajectory, PlotlyJS.scatter(x=[column2[1]], y=[column3[1]], 
                    type="scatter1d", mode="markers", marker=markerCSV, showlegend=false))
    end

    layout = PlotlyJS.Layout(xaxis_title = nameColumn2, yaxis_title = nameColumn3, title=title)
    myPlot = PlotlyJS.plot(trajectory, layout)
    if toPlot
        display(myPlot)
    end


    PlotlyJS.savefig(myPlot, saveFileName)
end