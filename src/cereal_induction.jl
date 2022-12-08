using CSV
include("data_generation.jl")

n_reps = 7; #replicates per genotype
species = categorical(["Triticale", "Wheat", "Barley", "Oats"]);
genotypes = categorical([1,2,3]); #three genotypes per species
treatments = categorical(["Control", "MeJA", "Insect"]);
colnames = [:Species, :Induction, :Genotype, :uniqueID];

fake_data = dataframeinit(treatments1=treatments, treatments2=species, groups=genotypes, replicates=n_reps, dfnames=colnames)
generated_params = generatedata(df= fake_data, idcol=:uniqueID, nested=:Genotype, seed=4444)
CSV.write("/home/isaac/Coding/Stats/cereal_induction/data/fake_data/fake_data.csv",fake_data)
CSV.write("data/fake_data/generated_params.csv", generated_params)

