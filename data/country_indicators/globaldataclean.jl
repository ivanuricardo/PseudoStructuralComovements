using DrWatson
@quickactivate
using CSV, DataFrames, JLD2

function cleancountry(country::String, df::DataFrame; logtransform::Bool=false, difftransform::Bool=false)
    country_data = filter(row -> row.Location == country, df)
    vals = logtransform ? log.(country_data.Value) : country_data.Value
    diffvalues = difftransform ? diff(vals) : vals
    return diffvalues
end

countries = ["USA", "CAN", "DEU", "FRA", "GBR"]

interestrate = CSV.read(datadir("./country_indicators/ir.csv"), DataFrame)
gdp = CSV.read(datadir("./country_indicators/gdp.csv"), DataFrame)
production = CSV.read(datadir("./country_indicators/production.csv"), DataFrame)
cpi = CSV.read(datadir("./country_indicators/cpi.csv"), DataFrame)

matdata = fill(NaN, 4, 5, 116)
for (i, country) in enumerate(countries)
    matdata[1, i, :] = cleancountry(country, gdp, logtransform=true)[1:116]
end
for (i, country) in enumerate(countries)
    matdata[2, i, :] = cleancountry(country, production, logtransform=true)[1:116]
end
for (i, country) in enumerate(countries)
    matdata[3, i, :] = cleancountry(country, interestrate)[1:116]
end
for (i, country) in enumerate(countries)
    matdata[4, i, :] = cleancountry(country, cpi)[1:116]
end

save("./globaldata.jld2", Dict("matdata" => matdata))
