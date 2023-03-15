using CSV
using MLDataUtils: rescale!
using Plots
using StatsPlots
using MCMCChains
using GLM
using MixedModels
using TuringGLM

### Combine 2019-2021 insured acreages of grain crops
### Data sourced from the Canadian Grain Commission https://www.grainscanada.gc.ca/en/grain-research/statistics/varieties-by-acreage/

include("acreage_sums.jl");
mylist = [
    "data/acreages/varieties-2019-c-en.csv",
    "data/acreages/varieties-2020-c-en.csv",
    "data/acreages/varieties-2021-c-en.csv"
];
acreages3years = compiledf(files=mylist)

totalacres = sumacres(data = acreages3years[:,[:grain_code, :variety, :acres]], groupon = :variety, vals = :acres, idcol = :grain_code)
CSV.write("data/acreages/varieties_three_year_ac.csv", totalacres)


### Generate Fake Data to validate a Turing Model
begin
include("data_generation.jl")
n_reps = 7; #replicates per genotype
species = categorical(["Triticale", "Wheat", "Barley", "Oats"]); #Grain crops
genotypes = categorical([1,2,3]); #three genotypes per species
treatments = categorical(["Control", "MeJA", "Insect"]);
colnames = [:Induction, :Species, :Genotype, :uniqueID];

fake_data = dataframeinit(treatments1=treatments, treatments2=species, groups=genotypes, replicates=n_reps, dfnames=colnames)
generated_params = generatedata(df= fake_data, idcol=:uniqueID, nested=:Genotype, seed=4444)
CSV.write("data/fake_data/fake_data.csv",fake_data)
CSV.write("data/fake_data/generated_params.csv", generated_params)
end

### Prep data for model use
include("turing_model.jl")
analysis_data = CSV.read("data/fake_data/fake_data.csv", DataFrame)
analysis_data_spread = spreadvars(df=analysis_data, treat_types=[:Species,:Induction], interaction=false)
#To incorporate Control & Barley into the intercept
begin
    analysis_data_spread = analysis_data_spread[:, Not([:treat1_Control, :treat2_Barley])]
end

#To center the data run this
begin
    cols_list = names(analysis_data_spread)
    scale_vals = DataFrame()
    ads_centered = rescalecols(df=analysis_data_spread, collist=Symbol.(cols_list[6:length(cols_list)]), centers = scale_vals)
end
adsc_idx = stringcoltoint(df=ad_spread, stringcol=:Genotype, intcol=:idx)
y_vals = Float64.(adsc_idx.response)
preds = Matrix(adsc_idx[:, 6:10])
idx = Int.(adsc_idx.idx)
my_model = randomintercept_regression(y_vals, preds, idx)
num_chains = 4
chains = sample(my_model, NUTS(0.6), MCMCThreads(), 1_000, num_chains)
summarystats(chains) |> DataFrame |> println
plt = plot(chains)

function changechainnames(chain, df, range1, range2)
    _namedict = Dict()
    if length(range1) != length(range2)
        println("Range lengths do no match")
    else
        _oldnames = String.(names(chain)[range1])
        _newnames = String.(names(df[:,range2]))
        for i in 1:length(range1)
            _namedict[_oldnames[i]] = _newnames[i]
        end
        _newchain = replacenames(chain, _namedict)
    end
    return _newchain
end

newnamesch = changechainnames(chains, adsc_idx, 2:8, 6:10)
namedchdf = DataFrame(summarystats(newnamesch))
usedparams = CSV.read("data/fake_data/generated_params.csv", DataFrame)


### Model didn't return the parameters I was expecting
### Check to see if my data generation was wrong
include("data_generation.jl")
fake_data = dataframeinit(treatments1=treatments, treatments2=species, groups=genotypes, replicates=n_reps, dfnames=colnames)
fd_spread = spreadvars(df=fake_data, treat_types=[:Species,:Induction], interaction=false)
df_with_resp, sim_params = generatedata(fd_spread, 5:11, :Genotype, 4444)

dfr_spread = spreadvars(df=df_with_resp, treat_types=[:Species,:Induction], interaction=false)
adsc_idx = stringcoltoint(df=dfr_spread, stringcol=:Genotype, intcol=:idx)
y_vals = Float64.(adsc_idx.response)
preds = Matrix(adsc_idx[:, 6:12])
idx = Int.(adsc_idx.idx)
new_model = randomintercept_regression(y_vals, preds, idx)
regression_chains = sample(new_model, NUTS(0.5), MCMCThreads(), 1_000, num_chains)
summarystats(regression_chains)

fit(MixedModel, @formula(response ~ Induction + Species + (1|Genotype)), analysis_data)

fm = @formula(response ~ Induction + Species + (1|Genotype))
model = turing_model(fm, analysis_data)
glmchains = sample(model, NUTS(), MCMCThreads(), 1_000, num_chains)
summarystats(glmchains) |> DataFrame |> println

### Plotting
plotting_data = CSV.read("data/fake_data/fake_data.csv", DataFrame)

mp = @df plotting_data groupedboxplot(:Induction, :response, group = :Species)
plot!(mp, xlabel = "Induction Treatment", ylabel = "Silicon Content (%)")
png(mp, "induction_plot")
test_df = DataFrame(var1 = repeat(["A", "B", "C"], inner = 3, outer =3), var2 = repeat(["D", "E", "F"], inner = 9, outer = 1), yvals = rand(Normal(0,1), 27))
tplot = @df test_df groupedboxplot(:var1, :yvals, group = :var2)


####

# Running a the above models on my real data

####

analysis_data = CSV.read("data/real_data/cleaned_si_df.csv", DataFrame)
ad_spread = spreadvars(df=analysis_data, treat_types=[:Species,:Induction], interaction=false)
begin
    ad_spread = ad_spread[:, Not([:treat2_Control, :treat1_Barley])]
end

#To center the data run this
begin
    cols_list = names(ad_spread)
    scale_vals = DataFrame()
    ads_centered = rescalecols(df=ad_spread, collist=[:Si_ppm], centers = scale_vals)
end
adsc_idx = stringcoltoint(df=ad_spread, stringcol=:Genotype, intcol=:idx)
y_vals = Float64.(adsc_idx.Si_ppm)
preds = Matrix(adsc_idx[:, 8:13])
idx = Int.(adsc_idx.idx)
my_model = randomintercept_regression(y_vals, preds, idx)
num_chains = 4
chains = sample(my_model, NUTS(0.6), MCMCThreads(), 1_000, num_chains)
summarystats(chains) |> DataFrame |> println
plt = plot(chains)

centered_data = rescalecols(df=analysis_data, collist=[:Si_ppm], centers = scale_vals)
lm1 = fit(MixedModel, @formula(Si_ppm ~ Induction * Species + (1|Genotype)), centered_data)

fm = @formula(Si_ppm ~ Induction + Species + (1|Genotype))
model = turing_model(fm, analysis_data)
glmchains = sample(model, NUTS(), MCMCThreads(), 1_000, num_chains)
summarystats(glmchains) |> DataFrame |> println

mp = @df analysis_data groupedboxplot(:Induction, :Si_ppm, group = :Species)
plot!(mp, xlabel = "Induction Treatment", ylabel = "Silicon Content (%)")
png(mp, "images/induction_plot")
test_df = DataFrame(var1 = repeat(["A", "B", "C"], inner = 3, outer =3), var2 = repeat(["D", "E", "F"], inner = 9, outer = 1), yvals = rand(Normal(0,1), 27))
tplot = @df test_df groupedboxplot(:var1, :yvals, group = :var2)