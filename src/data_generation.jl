#module DataGeneration

using DataFrames;
using Random;
using Distributions;
using CategoricalArrays;

##Need to generalize this out so it is more useful
# potentially supply all groups in a single list then do some magic to iterate through them and find the number of combos
function dataframeinit(;treatments1::AbstractVector, treatments2::AbstractVector, groups::AbstractArray, replicates::Int, dfnames::Vector)
    coltreat1 = repeat(treatments1, inner = 1* (length(groups)*length(treatments2)*replicates), outer = 1)
    coltreat2 = repeat(treatments2, inner = 1*(length(groups)*replicates), outer = 1*length(treatments1))
    colgroup = repeat(groups, inner= (1*replicates), outer = 1* length(treatments2) * length(treatments1))
    colid = shuffle(collect(1:(length(treatments1)*length(treatments2)*length(groups) * replicates)))
    colgroupconcat = string.(coltreat2, colgroup)
    dfout = DataFrame([coltreat1, coltreat2, colgroupconcat, colid], dfnames)
    return dfout
end 

function generatedata(;df::DataFrame, idcol::Symbol, nested::Symbol, seed::Int)
    #This function returns the list of the parameters, and adds fake results to the inputted dataframe 
    Random.seed!(seed);
    dfnoid = select(df, Not([idcol, nested]))
    σ_dist = Normal(0,1)
    temp_list = []
    coeff_list = []
    resp_list = []
    for i in 1:ncol(dfnoid)
        slice = dfnoid[!,i] 
        append!(temp_list, levels(slice))
        n_params = length(levels(slice))-1 # one category will be the reference, so length -1
        params = [0.0]
        append!(params,rand(n_params).-0.5 * 10)
        append!(coeff_list, params) # generate random parameter estimates, first entry is always 0 as it will be the baseline
    end
    append!(temp_list, levels(df[!,nested])) #add in the various genotypes
    append!(coeff_list, rand(length(levels(df[:,nested])))) #add in jitter for the genotypes
    temp = DataFrame([temp_list, coeff_list], [:parameter, :coeff]) # combine parameter names and estimates
    for i in 1:nrow(dfnoid) # loop over rows
        y_resp = [] # empty list for holding parameters
        for j in 1:ncol(dfnoid) # loop over columns
            #match the value in the main df (eg. "triticale") with the corresponding generated parameter value
            #store all the parameter vals that apply to that row in the y_resp list
            coeff_match = temp[temp.parameter .== dfnoid[i,j], :coeff]
            append!(y_resp, coeff_match)
        end
        #sum the stored parameters to get the final "phenotype" and add some error
        fakeobservation = sum(y_resp) + rand(σ_dist) #phenotype + error
        append!(resp_list, fakeobservation) #append "observation" to the list of generated y values
    end
    df.response = resp_list
    return temp
end
#end