using Turing;
using Random;
using Distributions;
using Statistics;
using StatsBase;
using MLDataUtils: rescale!;
using LinearAlgebra: I
using DataFrames;
using CSV;
using Plots, StatsPlots;
#This module defines the models and helper fxns for
# the induction experiment on Canadian cereal crops.
# The models will be hierarchical bayesian models using Turing.jl

# Model to account for inter-species differences
# yᵢ ~ MvNormal(uᵢ, σ)
# μᵢ = αᵢ + β₁ * x₁ + β₂ * x₂
# where x₁ is the species id and x₂ is the treatment code
# the intercept should vary by species
# σ = truncated(Normal(0,5), 0, 20) 
# revisit σ prior for Si measurements
# αᵢ ~ Normal(ā,σₐ) for j=1:number of genotypes
# ā ~ Normal(0.75, 0.25)
# prior for the average Si content (%) across all the species in the study
# σₐ ~ Exponential(a number)
# prior for variance among genotypes of the grasses
# βⱼ ~ Normal(β, σᵦ) for j = 1:number of genotypesSS
# β̄ ~ Normal(0, 5)
# prior for distribution of slopes
# σᵦ ~ Exponential(1)

function spreadvars(;df::DataFrame, treat_types::Vector{Symbol}, interaction::Bool) 
  _df = df
  #spread data to wide format
  length(treat_types) > 2 ? println("Fxn can only handle two treatment types") :
  if interaction == true
    for tr1 in unique(_df[:, treat_types[1]])
      _df[:, "treat1_$tr1"] = ifelse.(_df[:, treat_types[1]] .== tr1, 1, 0)
      for tr2 in unique(_df[:, treat_types[2]])
        _df[:, "treat2_$tr2"] = ifelse.(_df[:,treat_types[2]] .== tr2, 1,0)
        #Create columns for interactions
        _df[:, "$(tr1)x$(tr2)"] = ifelse.((_df[:,treat_types[1]] .== tr1 .&& _df[:,treat_types[2]] .== tr2), 1,0 )
      end
    end
  else 
    for tr1 in unique(_df[:, treat_types[1]])
      _df[:, "treat1_$tr1"] = ifelse.(_df[:, treat_types[1]] .== tr1, 1, 0)
      for tr2 in unique(_df[:, treat_types[2]])
        _df[:, "treat2_$tr2"] = ifelse.(_df[:,treat_types[2]] .== tr2, 1,0)
      end
    end
  end
  return _df
end

function rescalecols(;df::DataFrame, collist::Vector{Symbol}, centers::DataFrame)
  _df = df
  mus = []
  sigmas = []
  symbols = []
  centered_df_names = [:variable, :μ, :σ] 
  for i in collist
    _df[!,i] = convert.(Float64, _df[:, i])
    mu, sigma = rescale!(_df[!, i]) # rescale each predictor
    push!(mus, mu)
    push!(sigmas, sigma)
    push!(symbols, i)
  end
  append!(centers, DataFrame([symbols, mus, sigmas], centered_df_names))
  return _df
end

function stringcoltoint(;df::DataFrame, stringcol::Symbol, intcol::Symbol)
  _df = df
  _dict = Dict()
  _vec = Vector()
  #fill dict
  for (index, entry) in enumerate(unique(_df[:, stringcol]))
      _dict[entry] = index
  end
  #fill vec
  for strings in _df[:, stringcol]
      push!(_vec, _dict[strings])
  end
  _df[:, intcol] = _vec
  return _df
end
    
@model function randomintercept_regression(y,X,idx; n_gr=length(unique(idx)), predictors=size(X, 2))
  #priors, depends on data that has been rescaled
  α ~ Normal(0, 10) # population-level intercept
  β ~ filldist(Normal(0,10), predictors)
  σ ~ Exponential(10) # residual SD
  #prior for variance of random intercepts
  #usually requires thoughtful specification
  τ ~ truncated(Cauchy(0, 10); lower=0) #group level SDs
  αⱼ ~ filldist(Normal(0, τ), n_gr) # group level intercepts
  #likelihood
  ŷ = α .+ X * β .+ αⱼ[idx] 
  y ~ MvNormal(ŷ, σ^2 * I)
end