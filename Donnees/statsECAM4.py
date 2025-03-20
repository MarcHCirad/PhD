import pandas as pd

data = pd.read_stata("/home/hetier/Documents/StatCameroun/ecam4ind.dta")

print(data.shape)
dataSud = data[data["S0Q1"] == "Sud"]
dataSud = data.apply(lambda l: l["S0Q1"] == "Sud",axis=1)
print(dataSud.shape)
columns = data.columns
myColumns = ["S0Q1", "S0Q9", "S0Q11", "S05Q5", "S05Q20", "S10Q17", "S12Q1", "S12Q2", "S12Q3", "S12Q10"]
myColumns = myColumns + ["S02Q3A", "S02Q3B", "S02Q3C", "S12Q14", "S12Q54", "S12Q64"] + ["S12Q"+f"{k}" for k in range(71,94)]
selected_columns = []

for col_target in myColumns:
        if col_target in columns:
            selected_columns.append(col_target)
        else:
            print(col_target)
dataSudColumns = dataSud[selected_columns]

dataSudColumns.to_csv("/home/hetier/Documents/StatCameroun/ecam4indsort.csv", index=False)
