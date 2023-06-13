using CSV
using DataFrames
using Statistics
begin
#I first cut out each table of the plate data, see the attached plate1_matrix.csv file to see the format
#import each csv as a matrix, without the header
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

#make a list of the matrices (in R it should be c(p1m, p2m, ...))
matlist = [p1m, p2m, p3m, p4m, p5m, p6m, p7m, p8m, p9m, p10m]
#These are matrices of which samples were in each well. See attached file
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
#make another list
keylist = [p1k, p2k, p3k, p4k, p5k, p6k, p7k, p8k, p9k, p10k]
end
```
This function compiles a dataframe from a list of absorbance matrices and a list of key matrices

```


function builddf(ml, kl) 
    #make an empty list for keys
    keys = String[];
    #make an empty list for the absorbance values
    vals = Float64[];
    #error check to ensure an equal number of key and value matrices
    if length(ml) != length(kl)
        print("Unequal list lengths")
    else 
        #for each matrix is the list of value matrices        
        for m in eachindex(ml)
            #temporary matrix is the matrix at entry m in the value matrix list
            tempmat = ml[m]
            #temporary key matrix is the matrix at entry m in the key matrix list
            tempkey = kl[m]
            if size(tempkey) != size(tempmat)
                print("Error: Key and matrix different sizes")
            else
                #for each row
                for row in 1:(size(tempkey)[1])
                    #for each column
                    for col in 1:(size(tempkey)[2])
                        #if there is a blank cell fill with "MT" for empty
                        if ismissing(tempkey[row,col])
                            #append this value to the keys list
                            push!(keys, "MT")
                            #append the value to the vals list
                            push!(vals, tempmat[row,col])
                        else
                            #append the key to the keylist
                            push!(keys, string(tempkey[row,col]))\
                            #append the value to the value list
                            push!(vals, tempmat[row,col])

                            #print(typeof(keys))
                            #print(typeof(tempkey[row,col]))
                        end
                    end
                end
            end
        end
    end
    #create a dataframe of the keys list and the vals list
    return DataFrame(uniqueID = keys, absorbance = vals)
end

#now actually do the function above
absorbances = builddf(matlist, keylist)
#drop the MT keys
ab1 = filter(:uniqueID => k -> k != "MT", absorbances)
#drop the blanks
ab2 = filter(:uniqueID => k -> k != "BL", ab1)
#drop the gallic acid cells
ab3 = filter(:uniqueID => k -> k != "GA", ab2)
#drop the one we didn't know what it was
ab4 = filter(:uniqueID => k -> k != "?", ab3)
#group by uniqueID
gdf = groupby(ab4, :uniqueID)
#average at the uniqueID level to generate a single absorbance value per sample
#also count the number of rows used for the average
summ = combine(gdf, :absorbance => mean => :mean_absorbance, nrow)
#i had some doubles of samples, so if the number of rows > 6 I excluded it
final_df = filter(:nrow => n -> n!= 6, summ)
#save the absorbances dataframe
CSV.write("data/real_data/absorbances.csv", final_df)
#read it back in
finaldf = CSV.read("data/real_data/absorbances.csv", DataFrame) 
#read in the spreadsheet of powder masses                       
massdf = CSV.read("data/real_data/powdermass.csv", DataFrame)
#combine them using uniqueID as the key
finalmassdf = innerjoin(finaldf, massdf, on = :uniqueID)
#make a mass-corrected absorbance value
finalmassdf.mcabsorbance = finalmassdf.mean_absorbance ./ finalmassdf.powdermassg
finalmassdf
#save the df
CSV.write("data/real_data/mass_corrected_absorbances.csv", finalmassdf)

            

