module TuringModel

using #imports

#This module defines the models and helper fxns for
# the induction experiment on Canadian cereal crops.
# The models will be hierarchical bayesian models using Turing.jl

# Model to account for inter-species differences
# yᵢ ~ MvNormal(uᵢ, σ)
# μᵢ = αᵢ + βᵢ * x
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

end