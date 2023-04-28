using CSV
using DataFrames
using Statistics
begin
p1m = CSV.File("data/real_data/plate1_matrix.csv", header = false) |> Tables.matrix
p2m = CSV.File("data/real_data/plate2_matrix.csv", header = false) |> Tables.matrix
p3m = CSV.File("data/real_data/plate3_matrix.csv", header = false) |> Tables.matrix
p4m = CSV.File("data/real_data/plate4_matrix.csv", header = false) |> Tables.matrix
p5m = CSV.File("data/real_data/plate5_matrix.csv", header = false) |> Tables.matrix
p6m = CSV.File("data/real_data/plate6_matrix.csv", header = false) |> Tables.matrix
p7m = CSV.File("data/real_data/plate7_matrix.csv", header = false) |> Tables.matrix
p8m = CSV.File("data/real_data/plate8_matrix.csv", header = false) |> Tables.matrix
p9m = CSV.File("data/real_data/plate9_matrix.csv", header = false) |> Tables.matrix
p10m = CSV.File("data/real_data/plate10_matrix.csv", header = false) |> Tables.matrix

matlist = [p1m, p2m, p3m, p4m, p5m, p6m, p7m, p8m, p9m, p10m]

p1k = CSV.File("data/real_data/plate1_key.csv", header = false) |> Tables.matrix
p2k = CSV.File("data/real_data/plate2_key.csv", header = false) |> Tables.matrix 
p3k = CSV.File("data/real_data/plate3_key.csv", header = false) |> Tables.matrix
p4k = CSV.File("data/real_data/plate4_key.csv", header = false) |> Tables.matrix
p5k = CSV.File("data/real_data/plate5_key.csv", header = false) |> Tables.matrix
p6k = CSV.File("data/real_data/plate6_key.csv", header = false) |> Tables.matrix
p7k = CSV.File("data/real_data/plate7_key.csv", header = false) |> Tables.matrix
p8k = CSV.File("data/real_data/plate8_key.csv", header = false) |> Tables.matrix
p9k = CSV.File("data/real_data/plate9_key.csv", header = false) |> Tables.matrix
p10k = CSV.File("data/real_data/plate10_key.csv", header = false) |> Tables.matrix

keylist = [p1k, p2k, p3k, p4k, p5k, p6k, p7k, p8k, p9k, p10k]
end
function builddf(ml, kl)
    keys = String[];
    vals = Float64[];
    if length(ml) != length(kl)
        print("Unequal list lengths")
    else 
        for m in eachindex(ml)
            tempmat = ml[m]
            tempkey = kl[m]
            if size(tempkey) != size(tempmat)
                print("Error: Key and matrix different sizes")
            else
                for row in 1:(size(tempkey)[1])
                    for col in 1:(size(tempkey)[2])
                        if ismissing(tempkey[row,col])
                            push!(keys, "MT")
                            push!(vals, tempmat[row,col])
                        else
                            push!(keys, string(tempkey[row,col]))
                            push!(vals, tempmat[row,col])

                            #print(typeof(keys))
                            #print(typeof(tempkey[row,col]))
                        end
                    end
                end
            end
        end
    end
    return DataFrame(uniqueID = keys, absorbance = vals)
end

absorbances = builddf(matlist, keylist)
ab1 = filter(:uniqueID => k -> k != "MT", absorbances)
ab2 = filter(:uniqueID => k -> k != "BL", ab1)
ab3 = filter(:uniqueID => k -> k != "GA", ab2)
ab4 = filter(:uniqueID => k -> k != "?", ab3)
gdf = groupby(ab4, :uniqueID)
summ = combine(gdf, :absorbance => mean => :mean_absorbance, nrow)
final_df = filter(:nrow => n -> n!= 6, summ)

CSV.write("data/real_data/absorbances.csv", final_df)
finaldf = CSV.read("data/real_data/absorbances.csv", DataFrame)                        
massdf = CSV.read("data/real_data/powdermass.csv", DataFrame)
finalmassdf = innerjoin(finaldf, massdf, on = :uniqueID)
finalmassdf.mcabsorbance = finalmassdf.mean_absorbance ./ finalmassdf.powdermassg
finalmassdf
CSV.write("data/real_data/mass_corrected_absorbances.csv", finalmassdf)

            

