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

function spreadvars(;df::DataFrame, treat_types::Vector{Symbol}, interaction::Bool) 
    _df = df
    #spread data to wide format
    length(treat_types) > 2 ? println("Fxn can only handle two treatment types") :
    for tr1 in unique(_df[:, treat_types[1]])
      _df[:, "treat1_$tr1"] = ifelse.(_df[:, treat_types[1]] .== tr1, 1, 0)
      for tr2 in unique(_df[:, treat_types[2]])
        _df[:, "treat2_$tr2"] = ifelse.(_df[:,treat_types[2]] .== tr2, 1,0)
        #Create columns for interactions
        if interaction
          _df[:, "$(tr1)x$(tr2)"] = ifelse.((_df[:,treat_types[1]] .== tr1 .&& _df[:,treat_types[2]] .== tr2), 1,0 )
        end
      end
    end
    return _df
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
    println(length(resp_list))
    println(nrow(temp))
    df.response = resp_list
    return temp
end

function generatedata(df::DataFrame, range::UnitRange{Int64}, groups::Symbol, seed::Int)
    Random.seed!(seed)
    #y = α x βX + αⱼ + ϵ
    #αⱼ = intercept matrix * intercept values 
    _df = df
    gr_range = [] # saves the names of the intercept parameters
    for gr in unique(_df[:, groups])
        df[:, "$gr"] = ifelse.(_df[:, groups] .== gr, 1, 0)
        push!(gr_range, Symbol("$gr"))
    end
    intercepts = select(df, gr_range) #intercept matrix
    α_dist = Normal(0,1)
    α_list = rand(α_dist, size(intercepts, 2)) #intercept values
    preds = _df[:, range] #aka X
    params_list = names(preds)
    intercept_names = names(intercepts)

    σ_dist = Normal(0,1)
    σ_list = rand(σ_dist, size(preds, 1)) #aka ϵ
    coeff_dist = Normal(0,3)
    coeff_list = rand(coeff_dist, size(preds, 2)) * 1.25 # aka β
    α = rand(1)*5
    _df.response = Matrix(preds) * coeff_list + σ_list + (Matrix(intercepts) * α_list) .+ α[1]
    pushfirst!(params_list, "α")
    pushfirst!(coeff_list, α[1])
    append!(params_list, intercept_names)
    append!(coeff_list, α_list)
    paramsdf = DataFrame([params_list, coeff_list], [:parameter, :coefficient])
    select!(_df, Not(range))
    select!(_df, Not(gr_range))
    return _df, paramsdf
end


#end