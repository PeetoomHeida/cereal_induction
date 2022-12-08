#module DataGeneration

using DataFrames;
using CSV;
using Random;
using Distributions;
using CategoricalArrays;


n_reps = 7; #replicates per genotype
species = categorical(["Triticale", "Wheat", "Barley", "Oats"]);
genotypes = categorical([1,2,3]); #three genotypes per species
treatments = categorical(["Control", "MeJA", "Insect"]);
colnames = [:Species, :Induction, :Genotype, :uniqueID];
#should be 7 * 4 * 3 * 3 rows = 252 rows
function dataframeinit(;treatments1::AbstractVector, treatments2::AbstractVector, groups::AbstractArray, replicates::Int, dfnames::Vector, seed::Int)
    Random.seed!(seed)
    coltreat1 = repeat(treatments1, inner = 1* (length(groups)*length(treatments2)*replicates), outer = 1)
    coltreat2 = repeat(treatments2, inner = 1*(length(groups)*replicates), outer = 1*length(treatments1))
    colgroup = repeat(groups, inner= (1*replicates), outer = 1* length(treatments2) * length(treatments1))
    colid = shuffle(collect(1:(length(treatments1)*length(treatments2)*length(groups) * replicates)))
    dfout = DataFrame([coltreat1, coltreat2, colgroup, colid], dfnames)
    return dfout
end 


temp = dataframeinit(treatments1=species, treatments2=treatments, groups=genotypes, replicates=n_reps, dfnames=names, seed=1996)
m = @which dataframeinit(species, treatments, genotypes, n_reps, names)
Base.delete_method(m)

testdf = DataFrame([categorical(['a','b','c']), categorical(['d','e','f']), [1,2,3]], [:spp, :trt, :uniqueID])
outdf = DataFrame(colnames = Any[])
function generatedata(;df::DataFrame, idindex::Symbol, seed::Int, ydist::Distribution)
    Random.seed!(seed);
    dfnoid = select(df, Not(idindex))
    
    #Loop over each row in the dataframe
    temp = DataFrame(colname = Any[], levels = Any[], nlevels = [], coeffs = Any[])
    for i in ncol(dfnoid)
        slice = dfnoid[!,i] 
        push!(temp, [
        names(dfnoid)[i], #column name
        levels(slice), #levels in the column
        length(levels(slice)), # number of levels
        rand(length(levels(slice))),
         ]) # coefficients for the columns
    end
    return temp
    #= for i in nrow(df)
        for j in ncol(dfnoid)
            levels(dfnoid[!,j])
            temp
         

    end =#
end

generatedata(df= testdf, idindex=:uniqueID, seed=4444, ydist=Normal())


##end