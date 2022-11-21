module TuringModel
using Turing
using imports

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
# βⱼ ~ Normal(β, σᵦ) for j = 1:number of genotypes
# β̄ ~ Normal(0, 5)
# prior for distribution of slopes
# σᵦ ~ Exponential(1)

@model function multislope_regression(x, y, groups)
  
    n_gr = length(levels(groups)) # groups refers to cultivars, each cultivar is unique accross spp.
    #priors
    α ~ truncated(Normal(1, 0.75), 0, 15) # population-level intercept
    σ ~ Exponential(1) # residual SD
    #prior for variance of random intercepts
    #usually requires thoughtful specification
    τ ~ truncated(Normal(0, 5), 0, Inf)     # species-level SDs of intercepts
    αⱼ ~ filldist(Normal(0, τ), n_gr)       # species-level intercepts
    ζ ~ truncated(Normal(0,5), 0, Inf)    #species-level SDs slopes
    βⱼ ~ filldist(Normal(0,ζ), n_gr)  # species-level coefficients
    #likelihood
    ŷ = α .+ x .* βⱼ[groups] .+ αⱼ[groups] 
    y ~ MvNormal(ŷ, σ)
  end

end