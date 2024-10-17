import numpy as np

### Source : DESA RAPPORT CFSAM 2019, région Est 2019

#### Culture ###
productionsT = np.array([4840,18489,30232,840,892657,459043,841323])
superficieHa = np.array([400,19165,20249,105,79175,54343,89106])

rendementsTHa = productionsT / superficieHa
rendementskgHa = 1000 * rendementsTHa

proportionsCultivées = superficieHa / sum(superficieHa)

rendementskgHaMoyen = sum(proportionsCultivées * rendementskgHa)

print("Rendement moyen en kg par hectare : ", rendementskgHaMoyen)

#### Elevage ###
productionsViandeT = np.array([6747,492,1206,623,2531])

productionsViandeTHa = sum(productionsViandeT) / sum(superficieHa)
print("Production viande moyen en kg par hectare : ", productionsViandeTHa)