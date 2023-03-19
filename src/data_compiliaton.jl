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
xrf_data = CSV.read("data/real_data/xrf_data_86_cleaned.csv", DataFrame)
biomass_data = CSV.read("data/real_data/biomass_no106.csv", DataFrame)
xrf_colnames = names(xrf_data)
xrf_si_slice = xrf_data[:, ["Si Concentration", "Mn Concentration", "Sample ID"]] #select only the Si, Mn, and ID columns
rename!(xrf_si_slice, "Si Concentration" => "Si_ppm", "Sample ID" => "uniqueID", "Mn Concentration" => "Mn_ppm") #rename so you don't have to deal with spaces
combined_data_1 = innerjoin(herbivory_data, xrf_si_slice, on = :uniqueID) #inner join to move Si into the main dataframe
# The biomass data had two samples for 106, as did the XRF data. Unable to tell which was which so dropp it from the data set.
filter!(row -> row.uniqueID != 106, combined_data_1)
combined_data_2 = innerjoin(combined_data_1, biomass_data, on = :uniqueID)

function damageCat(x,y) 
    #This function lumps the Moderate and Severe Damage Classes into one "Damaged" category
    #and assigns "Undamaged" to all other plants in the "Insect" induction treatment
    #and assigns <missing> to all other plants
    if (x == "Insect" && y == "N")
        return "Undamaged"
    elseif (x == "Insect" && (y == "S" || y == "M"))
        return "Damaged"
    else
        return missing
    end
end

combined_data_2.isDamaged = damageCat.(combined_data_2.Induction, combined_data_2.Damage_rank)

#filtered_combined_data = filter(row -> (ismissing(row.Damage_rank) || row.Damage_rank != "N") && !ismissing(row.Si_ppm), combined_data )

#CSV.write("data/real_data/cleaned_si_df.csv", filtered_combined_data)
CSV.write("data/real_data/biomass_si_data.csv", combined_data_2)




