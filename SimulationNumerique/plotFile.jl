function plotTrajectory(inputFileName::String, saveFileName::String;
        plotType = "1d",
        variables = [],
        title = "",
        toPlot = false)
    myCSV = CSV.read(inputFileName, DataFrame)

    header = names(myCSV)
    time = myCSV[1:end, 1]
    nameTime = header[1]

    column2 = myCSV[2:end, 2]
    nameColumn2 = header[2]

    column3 = myCSV[2:end, 3]
    nameColumn3 = header[3]

    column4 = myCSV[2:end, 4]
    nameColumn4 = header[4]

    if plotType=="1d"
        l1 = PlotlyJS.Layout(yaxis_title = nameColumn2)
        s1 = PlotlyJS.scatter(x=time, y = column2)
        p1 = PlotlyJS.plot(s1, l1)

        l2 = PlotlyJS.Layout(yaxis_title = nameColumn3)
        s2 = PlotlyJS.scatter(x=time, y = column3)
        p2 = PlotlyJS.plot(s2, l2)

        l3 = PlotlyJS.Layout(xaxis_title = nameTime, yaxis_title = nameColumn4)
        s3 = PlotlyJS.scatter(x=time, y = column4)
        p3 = PlotlyJS.plot(s3, l3)

        myPlot = [p1 ; p2 ; p3]
        relayout!(myPlot, showlegend=false)
        if toPlot
            display(myPlot)
        end
        
    elseif plotType=="3d"
        l3d = PlotlyJS.Layout(scene = attr(xaxis_title = nameColumn2, yaxis_title=nameColumn3, zaxis_title=nameColumn4), title=title)
        s3d = PlotlyJS.scatter(x=column2, y=column3, z=column4, type="scatter3d", mode="lines",line=attr(color="darkblue", width=1))
        myPlot = PlotlyJS.plot(s3d, l3d)

        if toPlot
            display(myPlot)
        end
    end


    PlotlyJS.savefig(myPlot, saveFileName)
end

function plotBifurcationFile(inputFileName::String, saveFileName::String;
        xlabel = "",
        ylabel = "",
        title = "",    
        toPlot=false)
    myCSV = CSV.read(inputFileName, DataFrame)
    
    listParamX = tryparse.(Float64, myCSV[2:end, 1])
    listParamY = tryparse.(Float64, Vector(myCSV[1, 2:end]))

    stabilityValues = Matrix(myCSV[2:end, 2:end])

    dicStabilityNbr = Dict([("FVH",0), ("VH",1), ("FH",2), ("VH & FH",3), ("_", 4)])
    color = similar(stabilityValues, Int8)

    for ind_col in 1:length(listParamY)
        for ind_row in 1:length(listParamX)
            color[ind_col, ind_row] = dicStabilityNbr[stabilityValues[ind_row, ind_col]]
        end
    end
    
    layout = PlotlyJS.Layout(xaxis_title= xlabel, yaxis_title=ylabel, title=title)
    data = PlotlyJS.heatmap(x=listParamX, y=listParamY, z=color,
                            colorbar=attr(tickmode="array",
                                            tickvals=[0,1,2,3],
                                            ticktext=["EE^{FVH}", "EE^{VH}","EE^{FH}","EE^{FH} & EE^{VH}"]),
                                            # ticktext=[L"$EE^{FVH}$", L"$EE^{VH}$", L"$EE^{FH}$", L"$EE^{VH} & EE^{FH}$"]),
                            autocolorscale = false,
                            colorscale=[(0, "red"), (0.25, "red"), (0.25, "green"), (0.5, "green"), (0.5, "blue"), 
                                            (0.75, "blue"), (0.75, "black"), (0.99, "black"), (0.99, "white"), (1, "white")],
                            zmin = -0.5,
                            zmax = 3.5)

    myPlot = PlotlyJS.plot(data, layout)

    PlotlyJS.savefig(myPlot, saveFileName)

    if toPlot
        display(myPlot)
    end

end
