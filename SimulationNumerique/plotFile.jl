
function plotBifurcationFile(fileName::String)
    myCSV = CSV.read(fileName, DataFrame)
    listLambdaFH = tryparse.(Float64, myCSV[2:end, 1])
    listLambdaVH = tryparse.(Float64, Vector(myCSV[1, 2:end]))

    stabilityValues = Matrix(myCSV[2:end, 2:end])

    dicStabilityNbr = Dict([("FVH",0), ("VH",1), ("FH",2), ("VH & FH",3)])
    color = similar(stabilityValues, Int8)

    for ind_col in 1:length(listLambdaVH)
        for ind_row in 1:length(listLambdaFH)
            color[ind_row, ind_col] = dicStabilityNbr[stabilityValues[ind_row, ind_col]]
        end
    end
    
    layout = PlotlyJS.Layout(xaxis_title=L"$\lambda_{FH}$", yaxis_title=L"$\lambda_{VH}$")
    data = PlotlyJS.heatmap(x=listLambdaFH, y = listLambdaVH, z=color,
                            colorbar=attr(tickmode="array",
                                            tickvals=[0,1,2,3],
                                            ticktext=["EE^{FVH}", "EE^{VH}","EE^{FH}","EE^{FH} & EE^{VH}"]),
                                            # ticktext=[L"$EE^{FVH}$", L"$EE^{VH}$", L"$EE^{FH}$", L"$EE^{VH} & EE^{FH}$"]),
                            autocolorscale = false,
                            colorscale=[(0, "red"), (0.25, "red"), (0.25, "blue"), (0.5, "blue"), (0.5, "green"), 
                                            (0.75, "green"), (0.75, "black"), (1, "black")],
                            zmin = -0.5,
                            zmax = 3.5)
    plot = PlotlyJS.plot(data, layout)
    # display(plot)
    PlotlyJS.savefig(plot, "test.html")

end
