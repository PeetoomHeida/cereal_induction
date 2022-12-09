using CSV
using DataFrames

function compiledf(;files::AbstractVector, colnames::Vector{Symbol}) 
    temp = DataFrame()
    for i in 1:length(files)
        df = CSV.read(files[i], DataFrame)
        append!(temp, df)
    end
    rename!(temp, colnames)
    return temp
end;

function compiledf(;files::Vector{String}) 
    temp = DataFrame()
    for i in 1:length(files)
        df = CSV.read(files[i], DataFrame)
        append!(temp, df, promote=true)
    end
    return temp
end;

function sumacres(;data::DataFrame, groupon::Symbol, vals::Symbol, idcol::Symbol)
    internaldf = groupby(data, groupon) #group
    summeddf = combine(internaldf, vals => sum => vals) #sum
    sdfnames = leftjoin(summeddf, unique(data[:, [idcol, groupon]], groupon), on = groupon) #join idcol
    sort!(sdfnames, [idcol,vals], rev=true) #sort descending
end;