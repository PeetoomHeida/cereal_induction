using DataFrames
using GLM
using CSV

calibdata = CSV.read("./data/real_data/absorbances_calibration.csv", DataFrame)

fm = @formula (ppm ~ absorbance)
myreg = lm(fm, calibdata)
plot(calibdata.absorbance, calibdata.ppm, seriestype = :scatter)