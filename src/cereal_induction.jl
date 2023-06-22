begin
using CSV
using MLDataUtils: rescale!
using Plots
using StatsPlots
using MCMCChains
using GLM
using MixedModels
using TuringGLM
using DataFrames
using Measurements
using Statistics
using ColorSchemes
using Measures

function standarderror(col) 
    sd = std(col)
    denom = sqrt(length(col)-1)
    se = sd/denom
    return se
end
end
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

more_data = CSV.read("data/")
#############

# Running the above models on my real data

#############

include("linear_models.jl")
include("turing_model.jl")
analysis_data = CSV.read("data/real_data/biomass_si_data.csv", DataFrame)

#############

# Some summary stats on the data

#############
begin
my_names = unique(analysis_data.Genotype)
official_names = [
    "CS Camden", 
    "AC Summit", 
    "AC Morgan", 
    "CDC Austenson",
    "AAC Brandon",
    "CDC Landmark",
    "AB Stampeder",
    "Pronghorn",
    "CDC Plentiful",
    "CDC Copeland",
    "AAC Synergy",
    "Tyndal" ]

cultivar_dict = Dict(my_names .=> official_names)
end
begin
transform!(analysis_data, :Genotype => ByRow(x->cultivar_dict[x]) => :OfficialName)
analysis_data_noinduction = filter(row -> row.Induction == "Control", analysis_data)
analysis_data_noinduction.Si_ppm = analysis_data_noinduction.Si_ppm./10000
gdat_geno = groupby(analysis_data_noinduction, :Genotype)
summ_stats_geno = combine(gdat_geno, :Si_ppm => mean, :Si_ppm => std, nrow)
summ_stats_geno.Si_ppm_se = summ_stats_geno.Si_ppm_std ./ sqrt.(summ_stats_geno.nrow)
#CSV.write("./data/real_data/summary_statistics_genotype.csv", summ_stats_geno)
end
begin
gsi_spp = groupby(analysis_data_noinduction, [:OfficialName])
spp_si_contents_noind = combine(gsi_spp, [:Si_ppm, :Species] => ((si, sp) -> (sppMEANsi = mean(si), sppSEsi = standarderror(si), spp = first(sp))) => AsTable)
spp_si_contents_noind.GenotypeID = 1:12
sort!(spp_si_contents_noind, [:spp, :OfficialName])
end
spp_palette = palette([
    :palegreen, 
    :lime, 
    :darkgreen,
    :skyblue, 
    :dodgerblue, 
    :blue, 
    :pink, 
    :hotpink, 
    :deeppink,
    :orange, 
    :darkorange2, 
    :orange4])
begin
color_map = DataFrame(Cultivar = ["austenson", "brandon", "camden", "copeland", "landmark", "morgan", "plentiful", "pronghorn", "stampeder", "summit", "synergy", "tyndal"],
Order = [11, 4, 1, 10, 3, 6, 9, 8, 12, 5, 2, 7,],
Color = [
    :palegreen, 
    :lime, 
    :darkgreen,
    :skyblue, 
    :dodgerblue, 
    :blue, 
    :pink, 
    :hotpink, 
    :deeppink,
    :orange, 
    :darkorange2, 
    :orange4])

sort!(color_map, :Order)
spp_palette = color_map.Color
end

begin
spp_mean_si_plot = plot(spp_si_contents_noind.OfficialName[1:3], 
spp_si_contents_noind.sppMEANsi[1:3] .± spp_si_contents_noind.sppSEsi[1:3],
seriestype = :scatter,
markersize = 10,
#legend = spp_si_contents_noind.Genotype
labels = spp_si_contents_noind.spp[1])
plot!(spp_si_contents_noind.OfficialName[4:6], 
spp_si_contents_noind.sppMEANsi[4:6] .± spp_si_contents_noind.sppSEsi[4:6],
seriestype = :scatter,
markersize = 10,
#legend = spp_si_contents_noind.Genotype
labels = spp_si_contents_noind.spp[4])
plot!(spp_si_contents_noind.OfficialName[7:9], 
spp_si_contents_noind.sppMEANsi[7:9] .± spp_si_contents_noind.sppSEsi[7:9],
seriestype = :scatter,
markersize = 10,
#legend = spp_si_contents_noind.Genotype
labels = spp_si_contents_noind.spp[7])
plot!(spp_si_contents_noind.OfficialName[10:12], 
spp_si_contents_noind.sppMEANsi[10:12] .± spp_si_contents_noind.sppSEsi[10:12],
seriestype = :scatter,
markersize = 10,
#legend = spp_si_contents_noind.Genotype
labels = spp_si_contents_noind.spp[10])
plot!(xrotation =45)
plot!(xlabel = "Cultivar", ylabel = "Leaf Silicon Content (%)")
plot!(size = (600,600), dpi = 600)
plot!(xtickfontsize =14 ,xlabelfontsize = 16, ytickfontsize = 14, ylabelfontsize = 16, legendfontsize = 12)
plot!(grid=false)
plot!(left_margin=5mm, bottom_margin=2mm)
end
savefig(spp_mean_si_plot, "./manuscript/images/spp_si_content.png")


gdat_spp = groupby(analysis_data_noinduction, :Species)
summ_stats_spp = combine(gdat_spp, :Si_ppm => mean, :Si_ppm => std, nrow)
summ_stats_spp.Si_ppm_se = summ_stats_spp.Si_ppm_std ./ sqrt.(summ_stats_spp.nrow)
CSV.write("./data/real_data/summary_statistics_species.csv", summ_stats_spp)

si_mass_model = biomass_si_regression(df=analysis_data, interaction = false)
begin
    biomass_si_scatter = plot(analysis_data.mass_g, analysis_data.Si_ppm/10000, group = analysis_data.Species, seriestype=:scatter, smooth=true)
    plot!(biomass_si_scatter, xlabel = "Aboveground Biomass (g)", ylabel = "Leaf Si content (%)")
    plot!(size = (800,600))
    png(biomass_si_scatter, "images/biomass_regression")
end


scale_vals = DataFrame()
ads_centered = rescalecols(df=analysis_data, collist=[:Si_ppm], centers = scale_vals)
full_model = fullcenteredglm(df=analysis_data)
my_turing_model, center_vals = HMturingmodel(df = analysis_data, interaction = true)
insect_model = insectonlymodel(df=analysis_data, interaction=true)
insect_slice = analysis_data[(analysis_data.Induction .== "Insect"), :]
plot(insect_slice.isDamaged, insect_slice.Si_ppm, seriestype=:boxplot)
centered_data = rescalecols(df=analysis_data, collist=[:Si_ppm], centers = scale_vals)
lm1 = fit(MixedModel, @formula(Si_ppm ~ Induction * Species + mass_g + (1|Genotype)), centered_data)

fm = @formula(Si_ppm ~ Induction + Species + (1|Genotype))
model = turing_model(fm, analysis_data)
glmchains = sample(model, NUTS(), MCMCThreads(), 1_000, num_chains)
summarystats(glmchains) |> DataFrame |> println

#Plot [Si] against biomass, looks like a negative relationship
biomass_si_regression = lm(@formula(Si_ppm ~ mass_g * Species), analysis_data)
gboxplot = @df analysis_data groupedboxplot(:Induction, :Si_ppm/10000, group = :Species)
gboxplot = @df analysis_data groupeddotplot(:Induction, :Si_ppm/10000, group = :Species)
plot!(analysis_data.Induction, analysis_data.Si_ppm/10000, group = analysis_data.Species, seriestype=:scatter)
plot!(gboxplot, xlabel = "Induction Treatment", ylabel = "Silicon Content (%)", size = (800,600), dpi = 600)
png(gboxplot, "manuscript/images/induction_plot")
test_df = DataFrame(var1 = repeat(["A", "B", "C"], inner = 3, outer =3), var2 = repeat(["D", "E", "F"], inner = 9, outer = 1), yvals = rand(Normal(0,1), 27))
tplot = @df test_df groupedboxplot(:var1, :yvals, group = :var2)

filtered_analysis_data = filter(row -> ismissing(row.isDamaged) || row.isDamaged  ≠ "Undamaged", analysis_data)

g_analysisdata = groupby(filtered_analysis_data, [:Induction, :Species])
si_induction_summary = combine(g_analysisdata, [:Si_ppm, :Species, :Induction] => ((si, sp, in) -> (grpMEANsi = mean(si/10000), grpSEsi = standarderror(si/10000), spp = first(sp), ind = first(in))) => AsTable)


sort!(si_induction_summary, [:ind, :spp])
si_induction_summary.order = [1,2,3,4,5,6,7,8,9,10,11,12]
begin
    induction_si_plot = plot(si_induction_summary.order[[1,5,9]], 
    si_induction_summary.grpMEANsi[[1,5,9]] .± si_induction_summary.grpSEsi[[1,5,9]],
    seriestype = :scatter,
    markersize = 10,
    #legend = si_induction_summary.Genotype
    labels = "Barley")

    plot!(si_induction_summary.order[[2,6,10]], 
    si_induction_summary.grpMEANsi[[2,6,10]] .± si_induction_summary.grpSEsi[[2,6,10]],
    seriestype = :scatter,
    markersize = 10,
    #legend = si_induction_summary.Genotype
    labels = "Oats")

    plot!(si_induction_summary.order[[3,7,11]], 
    si_induction_summary.grpMEANsi[[3,7,11]] .± si_induction_summary.grpSEsi[[3,7,11]],
    seriestype = :scatter,
    markersize = 10,
    #legend = si_induction_summary.Genotype
    labels = "Triticale")

    plot!(si_induction_summary.order[[4,8,12]], 
    si_induction_summary.grpMEANsi[[4,8,12]] .± si_induction_summary.grpSEsi[[4,8,12]],
    seriestype = :scatter,
    markersize = 10,
    #legend = si_induction_summary.Genotype
    labels = "Wheat")

    #= plot!(si_induction_summary.spp[10:12], 
    si_induction_summary.grpMEANsi[10:12] .± si_induction_summary.grpSEsi[10:12],
    seriestype = :scatter,
    markersize = 8,
    #legend = si_induction_summary.Genotype
    labels = si_induction_summary.ind[10]) =#

    plot!(xticks = ([2.5, 6.5, 10.5], unique(si_induction_summary.ind)))
    plot!(xtickfontsize =16 ,xlabelfontsize = 20, ytickfontsize = 16, ylabelfontsize = 20, legendfontsize = 16)
    plot!(left_margin=5mm, bottom_margin=2mm)
    #plot!(xrotation =45)

    plot!(xlabel = "Treatment Type", ylabel = "Leaf Silicon Content (%)")
    plot!(grid=false)
    plot!(size = (800,600), dpi = 600)
    end

savefig(induction_si_plot, "./manuscript/images/induction_plot_bytrt.png")


ad_spread = spreadvars(df=analysis_data, treat_types=[:Species,:Induction], interaction=true)
    ad_spread = ad_spread[:, Not([:treat2_Control, :treat1_Barley])]
    scale_vals = DataFrame()
    ads_centered = rescalecols(df=ad_spread, collist=[:Si_ppm], centers = scale_vals)
    adsc_idx = stringcoltoint(df=ads_centered, stringcol=:Genotype, intcol=:idx)
    y_vals = Float64.(adsc_idx.Si_ppm)
    preds = Matrix(adsc_idx[:, 8:13])
    idx = Int.(adsc_idx.idx)
    my_model = randomintercept_regression(y_vals, preds, idx)
    num_chains = 4
    chains = sample(my_model, NUTS(0.6), MCMCThreads(), 1_000, num_chains)
    summarystats(chains) |> DataFrame |> println
    plot(chains)
#begin
phenolic_data = CSV.read("data/real_data/si_absorbance_data.csv", DataFrame)
#The regression formula to convert absorbance to ppm is:
#   ppm = -222.142 + 2073.3(absorbance)
function phenolic_conversion(x)
    #to convert from absorbance to ppm
    ppm = -222.142 + 2073.3(x)
    #to convert from ppm to mg
    #mg = ppm*L
    mg = ppm * 0.003

    #To make the value mg/kg use * 10,
    # to make the value %dry matter, use *10000 (mg/100mg)
   # y = ((x -0.05 + 222.142) * 0.003)/2073.3*100000
    return mg
end
begin
transform!(phenolic_data, :mean_absorbance => ByRow(x-> phenolic_conversion(x)) =>  :phenolicsmg)
phenolic_data.phenolicsmgpperg = phenolic_data.phenolicsmg ./ phenolic_data.powdermassg
phenolics_filtered = filter(row -> ismissing(row.isDamaged) || row.isDamaged  ≠ "Undamaged", phenolic_data)
phe_summary_stats_g  = groupby(phenolics_filtered, :Species)
summ_stats_phe = combine(phe_summary_stats_g, :phenolicsmgpperg => mean, :phenolicsmgpperg => standarderror)
gboxplot_ab = @df phenolics_filtered groupedboxplot(:Induction, :phenolicsmgpperg, group = :Species)
plot!(xlabel = "Treatment Group", ylabel = "Phenolic Content", dpi = 600, size = (800,600))
png(gboxplot_ab, "manuscript/images/absorbance_boxplots")
phenolic_data.Si_pc = phenolic_data.Si_ppm./10000
phenolic_data.treatmentcombo = phenolic_data.Induction .* " " .* phenolic_data.Species
ph_si_scatter = plot(phenolic_data.Si_pc, phenolic_data.phenolicsmgpperg, group=phenolic_data.Species, seriestype=:scatter, markersize =8)
plot!(ph_si_scatter, xlabel = "Leaf Silicon Content (%)", ylabel = "Leaf Phenolic Content (mg/g)")
plot!(size = (800,600), dpi = 800)
plot!(grid=false)
plot!(xtickfontsize =16 ,xlabelfontsize = 20, ytickfontsize = 16, ylabelfontsize = 20, legendfontsize = 16)
    plot!(left_margin=5mm, bottom_margin=2mm)
png(ph_si_scatter, "manuscript/images/phenolic_silicon_regression")
end
begin
g_phenolicdata = groupby(phenolics_filtered, [:Induction, :Species])
phe_induction_summary = combine(g_phenolicdata, [:phenolicsmgpperg, :Species, :Induction] => ((abs, sp, in) -> (grpMEANabs = mean(abs), grpSEabs = standarderror(abs), spp = first(sp), ind = first(in))) => AsTable)


sort!(phe_induction_summary, [:ind, :spp])
phe_induction_summary.order = [1,2,3,4,5,6,7,8,9,10,11,12]
end
begin
    induction_phe_plot = plot(phe_induction_summary.order[[1,5,9]], 
    phe_induction_summary.grpMEANabs[[1,5,9]] .± phe_induction_summary.grpSEabs[[1,5,9]],
    seriestype = :scatter,
    markersize = 10,
    #legend = phe_induction_summary.Genotype
    labels = "Barley")

    plot!(phe_induction_summary.order[[2,6,10]], 
    phe_induction_summary.grpMEANabs[[2,6,10]] .± phe_induction_summary.grpSEabs[[2,6,10]],
    seriestype = :scatter,
    markersize = 10,
    #legend = phe_induction_summary.Genotype
    labels = "Oats")

    plot!(phe_induction_summary.order[[3,7,11]], 
    phe_induction_summary.grpMEANabs[[3,7,11]] .± phe_induction_summary.grpSEabs[[3,7,11]],
    seriestype = :scatter,
    markersize = 10,
    #legend = phe_induction_summary.Genotype
    labels = "Triticale")

    plot!(phe_induction_summary.order[[4,8,12]], 
    phe_induction_summary.grpMEANabs[[4,8,12]] .± phe_induction_summary.grpSEabs[[4,8,12]],
    seriestype = :scatter,
    markersize = 10,
    #legend = phe_induction_summary.Genotype
    labels = "Wheat")

    #= plot!(phe_induction_summary.spp[10:12], 
    phe_induction_summary.grpMEANsi[10:12] .± phe_induction_summary.grpSEsi[10:12],
    seriestype = :scatter,
    markersize = 8,
    #legend = phe_induction_summary.Genotype
    labels = phe_induction_summary.ind[10]) =#

    plot!(xticks = ([2.5,6.5,10.5], unique(phe_induction_summary.ind)))

    #plot!(xrotation =45)
    plot!(xtickfontsize =18 ,xlabelfontsize = 20, ytickfontsize = 18, ylabelfontsize = 20, legendfontsize = 16)
    plot!(left_margin=5mm, bottom_margin=2mm)
    plot!(xlabel = "Treatment Type", ylabel = "Leaf Phenolic Content (mg/g)")
    plot!(size = (800,600), dpi = 600)
    plot!(grid=false)
    
    end

savefig(induction_phe_plot, "./manuscript/images/phenolic_induction_plot_bytrt.png")

begin
gph = groupby(phenolic_data, [:Species])
test = combine(gph, [:phenolicsmgpperg, :Si_ppm] => ((a, s) -> (sppMEANabsorbance = mean(a), sppSEabsorbance = standarderror(a),sppMEANsi = mean(s), sppSEsi = standarderror(s))) => AsTable)


meansphsiscatter = plot(test.sppMEANsi .± test.sppSEsi, test.sppMEANabsorbance .± test.sppSEabsorbance, seriestype = "scatter", groups = test.Species, markersize = 15, palette = :Dark2_3)
plot!(meansphsiscatter, xlabel = "Leaf Silicon Content (%)", ylabel = "Leaf Phenolic Content")
plot!(size = (800,600), dpi = 800)
#png(ph_si_scatter, "manuscript/images/phenolic_silicon_regression")

gph_ind = groupby(phenolic_data, [:Induction])
test_ind = combine(gph_ind, [:phenolicsmgpperg, :Si_pc] => ((a, s) -> (sppMEANabsorbance = mean(a), sppSEabsorbance = standarderror(a),sppMEANsi = mean(s), sppSEsi = standarderror(s))) => AsTable)

meansphsiscatter_induction = plot(test_ind.sppMEANsi .± test_ind.sppSEsi, test_ind.sppMEANabsorbance .± test_ind.sppSEabsorbance, seriestype = "scatter", groups = test_ind.Induction, markersize = 15)
plot!(meansphsiscatter_induction, xlabel = "Leaf Silicon Content (%)", ylabel = "Leaf Phenolic Content")
plot!(size = (800,600), dpi = 800)
#png(ph_si_scatter, "manuscript/images/phenolic_silicon_regression")

#= gph_trt = groupby(phenolics_filtered, [:treatmentcombo])
test_trt = combine(gph_trt, [:phenolicsmgpperg, :Si_pc, :Species, :Induction] => ((ab, si, sp, in) -> (sppMEANabsorbance = mean(ab), sppSEabsorbance = standarderror(ab),sppMEANsi = mean(si), sppSEsi = standarderror(si), spp = first(sp), ind = first(in))) => AsTable)
meansphsiscatter_induction = plot(test_trt.sppMEANsi .± test_trt.sppSEsi, 
    test_trt.sppMEANabsorbance .± test_trt.sppSEabsorbance, 
    seriestype = "scatter", 
    groups = test_trt.treatmentcombo, 
    markersize = 15
) =#

plot(test_trt.treatmentcombo, test_trt.sppMEANsi .± test_trt.sppSEsi, seriestype = :scatter, group = test_trt.spp, size = (1200,600), dpi = 600, markersize = 10)
plot!(xlabel = "Species Treatment Combination", ylabel = "Mean Silicon Content (%)")

plot(test_trt.treatmentcombo, test_trt.sppMEANabsorbance .± test_trt.sppSEabsorbance, seriestype = :scatter, group = test_trt.spp, size = (1200,600), dpi = 600, markersize = 10)
plot!(xlabel = "Species Treatment Combination", ylabel = "Mean Silicon Content (%)")
end
begin
transform!(phenolic_data, :Genotype => ByRow(x->cultivar_dict[x]) => :OfficialName)
phenolics_data_noinduction = filter(row -> row.Induction == "Control", phenolic_data)
gdat_geno = groupby(phenolics_data_noinduction, :Genotype)
summ_stats_geno = combine(gdat_geno, :phenolicsmgpperg => mean, :phenolicsmgpperg => standarderror, nrow)
#CSV.write("./data/real_data/summary_statistics_genotype.csv", summ_stats_geno)
end
begin
gab_spp = groupby(phenolics_data_noinduction, [:OfficialName])
spp_phe_contents_noind = combine(gab_spp, [:phenolicsmgpperg, :Species] => ((ab, sp) -> (sppMEANab = mean(ab), sppSEab = standarderror(ab), spp = first(sp))) => AsTable)
spp_phe_contents_noind.GenotypeID = 1:12
sort!(spp_phe_contents_noind, [:spp, :OfficialName])

spp_palette = palette([
    :palegreen, 
    :lime, 
    :darkgreen,
    :skyblue, 
    :dodgerblue, 
    :blue, 
    :pink, 
    :hotpink, 
    :deeppink,
    :orange, 
    :darkorange2, 
    :orange4])

color_map = DataFrame(Cultivar = ["austenson", "brandon", "camden", "copeland", "landmark", "morgan", "plentiful", "pronghorn", "stampeder", "summit", "synergy", "tyndal"],
Order = [11, 4, 1, 10, 3, 6, 9, 8, 12, 5, 2, 7,],
Color = [
    :palegreen, 
    :lime, 
    :darkgreen,
    :skyblue, 
    :dodgerblue, 
    :blue, 
    :pink, 
    :hotpink, 
    :deeppink,
    :orange, 
    :darkorange2, 
    :orange4])

sort!(color_map, :Order)
spp_palette = color_map.Color
plot(spp_phe_contents_noind.GenotypeID, 
spp_phe_contents_noind.sppMEANab .± spp_phe_contents_noind.sppSEab, 
seriestype = :scatter, 
group = spp_phe_contents_noind.OfficialName, 
size = (1200,600), 
dpi = 600, 
markersize = 10,
grid=false,
palette = spp_palette
)
plot!(xticks = ([2,5,8,11], unique(spp_phe_contents_noind.spp)))
plot!(xlabel = "Cultivar", ylabel = "Mean Silicon Content (%)")
end
begin
spp_mean_phe_plot = plot(spp_phe_contents_noind.OfficialName[1:3], 
spp_phe_contents_noind.sppMEANab[1:3] .± spp_phe_contents_noind.sppSEab[1:3],
seriestype = :scatter,
markersize = 10,
#legend = spp_phe_contents_noind.Genotype
labels = spp_phe_contents_noind.spp[1])
plot!(spp_phe_contents_noind.OfficialName[4:6], 
spp_phe_contents_noind.sppMEANab[4:6] .± spp_phe_contents_noind.sppSEab[4:6],
seriestype = :scatter,
markersize = 10,
#legend = spp_phe_contents_noind.Genotype
labels = spp_phe_contents_noind.spp[4])
plot!(spp_phe_contents_noind.OfficialName[7:9], 
spp_phe_contents_noind.sppMEANab[7:9] .± spp_phe_contents_noind.sppSEab[7:9],
seriestype = :scatter,
markersize = 10,
#legend = spp_phe_contents_noind.Genotype
labels = spp_phe_contents_noind.spp[7])
plot!(spp_phe_contents_noind.OfficialName[10:12], 
spp_phe_contents_noind.sppMEANab[10:12] .± spp_phe_contents_noind.sppSEab[10:12],
seriestype = :scatter,
markersize = 10,
#legend = spp_phe_contents_noind.Genotype
labels = spp_phe_contents_noind.spp[10])
plot!(xrotation =45)
plot!(xlabel = "Cultivar", ylabel = "Leaf Phenolic Content (mg/kg)")

plot!(size = (600,600), dpi = 600)
plot!(xtickfontsize =14 ,xlabelfontsize = 16, ytickfontsize = 14, ylabelfontsize = 16, legendfontsize = 12)
plot!(grid=false)
plot!(left_margin=5mm, bottom_margin=2mm)
end
savefig(spp_mean_phe_plot, "./manuscript/images/spp_phenolic_content.png")
