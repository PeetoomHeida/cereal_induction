using GLM
using MixedModels
using DataFrames


"""
This model is used to test for differences in the silicon content of plants in the insect treatment.

It filters the inputted DataFrame (in the same format as the CSV) to only select the Insect treated
plants.

Depending on the interaction parameter, it then runs one of these two 
linear mixed effects models: 

    Si_ppm ~ Species + Damage + mass_g + (1|Genotype) #interaction = false

    Si_ppm ~ Species * Damage + mass_g + (1|Genotype) #interaction = false

"""
function insectonlymodel(;df::DataFrame, interaction::Bool) 
    df1 = filter(row -> row.Induction == "Insect", df)
    if interaction 
        interactivemodel = fit(MixedModel, @formula(Si_ppm ~ Species * isDamaged + mass_g + (1|Genotype)), df1)
        return interactivemodel
    else
        additivemodel = fit(MixedModel, @formula(Si_ppm ~ Species + isDamaged + mass_g + (1|Genotype)), df1)
        return additivemodel
    end
end

"""
This model is used to test for an effect of our induction treatments on the silicon content of the plants.

It filters out the non-damaged plants in the insect treatment group, as these plants received a minimal amount of damage.

It then runs a model with the following formula:

    Si_ppm ~ Species * Induction + mass_g + (1|Genotype)
"""
function fullcenteredglm(;df::DataFrame)
    df1 = filter(row -> (ismissing(row.isDamaged)|| row.isDamaged != "Undamaged"), df)
    outputmodel = fit(MixedModel, @formula(Si_ppm ~ Species * Induction + mass_g + (1|Genotype)), df1)
    return outputmodel
end
"""
This model tests for a significat correlation between biomass and silicon content.
It can test the model with/without an interaction term between the species and biomass treatments

If interaction = true the following formula is used:

    Si_ppm ~ mass_g * Species + (1|Genotype))

If interaction = false the following formula is used:

    Si_ppm ~ mass_g + Species + (1|Genotype))
"""
function biomass_si_regression(;df::DataFrame, interaction::Bool)
    if interaction
        outputmodel = fit(MixedModel, @formula(Si_ppm ~ mass_g * Species + (1|Genotype)), df)
    else
        outputmodel = fit(MixedModel, @formula(Si_ppm ~ mass_g + Species + (1|Genotype)), df)
    end
    return outputmodel
end

"""
This function takes the data frame and preps it for a turing model.
Prep includes scaling/centering the predictor and response variables, and spreading the predictors into dummy coded variables.
It then calls the turing model and runs. The return value is the output of the sampling, an MCMCChains object.
It also returns the dataframe object for the scaled values of the predictors.
Be sure to assign the output to two objects. 
output: MCMCChains, DataFrame
"""
function HMturingmodel(;df::DataFrame, interaction::Bool)
    include("turing_model.jl")
    ad_spread = spreadvars(df=df, treat_types=[:Species,:Induction], interaction=interaction)
    ad_spread = ad_spread[:, Not([:treat2_Control, :treat1_Barley])]
    scale_vals = DataFrame()
    ads_centered = rescalecols(df=ad_spread, collist=[:Si_ppm], centers = scale_vals)
    adsc_idx = stringcoltoint(df=ads_centered, stringcol=:Genotype, intcol=:idx)
    y_vals = Float64.(adsc_idx.Si_ppm)
    if interaction
        preds = Matrix(adsc_idx[:, 12:28])
    else
        preds = Matrix(adsc_idx[:,12:16])
    end
    idx = Int.(adsc_idx.idx)
    my_model = randomintercept_regression(y_vals, preds, idx)
    num_chains = 4
    chains = sample(my_model, NUTS(0.6), MCMCThreads(), 1_000, num_chains)
    summarystats(chains) |> DataFrame |> println
    plot(chains)
    return chains, scale_vals
end

