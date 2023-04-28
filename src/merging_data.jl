using CSV
using DataFrames

abdata = CSV.read("data/real_data/mass_corrected_absorbances.csv", DataFrame)
sidata = CSV.read("data/real_data/biomass_si_data.csv", DataFrame)

merged = innerjoin(sidata, abdata, on=:uniqueID)

CSV.write("data/real_data/si_absorbance_data.csv", merged)