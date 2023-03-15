using CSV
using DataFrames

#=
This script is written to combine herbivory data with the raw output of the pXRF analysis. We generated herbivory data by observing
which cricket-exposed plants had suffered damage after the 18h exposure period. Damage was categorized into three classes: 
undamaged, moderately damaged (leaf chewed by not severed), severly damaged (leaf severed). The pXRF data were generated
using an Olympus Vanta pXRF device mounted in a benchtop stand. Prior to analysis, we pressed leaf tissue powder into a 13mm
pellet using a hydraulic press and die set. We pressed the powder to 300 bar of pressue, resulting in a stable pellet with a 
uniform surface. To quantify silicon, we used a 45 second scan time, a time that balanced throughput with accuracy (increasing 
scan time by a factor of two resulted in only a 7% reduction in StDev of the measurement).
=#

herbivory_data = CSV.read("data/real_data/induction_data.csv", DataFrame)
xrf_data = CSV.read("data/real_data/xrf_data.csv", DataFrame)
xrf_colnames = names(xrf_data)
xrf_si_slice = xrf_data[:, ["Si Concentration", "Sample ID"]] #select only the Si and ID columns
rename!(xrf_si_slice, "Si Concentration" => "Si_ppm", "Sample ID" => "uniqueID") #rename so you don't have to deal with spaces
combined_data = leftjoin(herbivory_data, xrf_si_slice, on = :uniqueID) #left join to move Si into the main dataframe
filtered_combined_data = filter(row -> (ismissing(row.Damage_rank) || row.Damage_rank != "N") && !ismissing(row.Si_ppm), combined_data )

CSV.write("data/real_data/cleaned_si_df.csv", filtered_combined_data)




