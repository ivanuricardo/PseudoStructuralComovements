using DrWatson
@quickactivate
using CSV, DataFrames, JLD2

function cleancountry(country::String, df::DataFrame; logtransform::Bool=false, diffval::Int=1)
    country_data = filter(row -> row.Location == country, df)
    vals = logtransform ? log.(country_data.Value) : country_data.Value
    if diffval > 0
        diffvalues = vals[diffval+1:end] .- vals[1:end-diffval]
    end
    return diffvalues
end

function seasonal_diff_log(x::Vector, s::Int=4)
    return log.(x[s+1:end]) .- log.(x[1:end-s])
end

countries = ["USA", "CAN", "DEU", "FRA", "GBR"]
#=countries = ["USA", "CAN", "DEU", "FRA"]=#

interestrate = CSV.read(datadir("./country_indicators/ir.csv"), DataFrame)
gdp = CSV.read(datadir("./country_indicators/gdp.csv"), DataFrame)
production = CSV.read(datadir("./country_indicators/production.csv"), DataFrame)
cpi = CSV.read(datadir("./country_indicators/cpi.csv"), DataFrame)

# Starts at 1995-Q1 until 2019-Q4, for 99 observations
matdata = fill(NaN, 4, 5, 99)
for (i, country) in enumerate(countries)
    matdata[1, i, :] = cleancountry(country, gdp, logtransform=true, diffval=1)[17:115]
end
for (i, country) in enumerate(countries)
    matdata[2, i, :] = cleancountry(country, production, logtransform=true, diffval=1)[17:115]
end
for (i, country) in enumerate(countries)
    matdata[3, i, :] = cleancountry(country, interestrate, diffval=1)[17:115]
end

for (i, country) in enumerate(countries)
    matdata[4, i, :] = cleancountry(country, cpi, logtransform=true, diffval=4)[14:112]
end

save("./globaldata.jld2", Dict("matdata" => matdata))
