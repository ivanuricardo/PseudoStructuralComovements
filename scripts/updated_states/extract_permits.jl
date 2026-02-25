using DrWatson
@quickactivate :PseudoStructuralComovements
using Downloads, CSV, DataFrames
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

base_url = "https://www2.census.gov/econ/bps/State/"
save_dir = datadir("updated_states/raw-permits/")
num_files = length(filter(f -> isfile(joinpath(save_dir, f)), readdir(save_dir)))
states = ["Iowa", "Illinois", "Indiana", "Minnesota", "North Dakota", "Ohio", "South Dakota"]

# for year in 2001:2019
#     yy = lpad(string(year % 100), 2, '0')
#     for month in 1:12
#         mm = lpad(string(month), 2, '0')
#         filename = "st$(yy)$(mm)c.txt"
#         url = base_url * filename
#         save_path = joinpath(save_dir, filename)
#         try
#             Downloads.download(url, save_path)
#             println("Downloaded $filename")
#         catch e
#             println("Failed to download $filename: $e")
#         end
#     end
# end

# 2015 feb doesn't seem to work, so i input it manually
# the yearly cumulative sum is available, so i use feb - jan for this value
# Remove the commented out version to reinstall the data
missing_dir = datadir("updated_states/missing_dir/")
# Downloads.download(base_url * "st1501y.txt", missing_dir * "st1501y.txt")
# Downloads.download(base_url * "st1502y.txt", missing_dir * "st1502y.txt")

jan = CSV.read(missing_dir * "st1501y.txt", DataFrame; header = 1)
feb = CSV.read(missing_dir * "st1502y.txt", DataFrame; header = 1)
missing_vec = fill(NaN, length(states))

for (i, state) in enumerate(states)
    jan_row = jan[(.!ismissing.(jan.State)) .& (jan.State .== state), :]
    feb_row = feb[(.!ismissing.(feb.State)) .& (feb.State .== state), :]

    one_unit_jan = parse(Float64, jan_row[1,7])
    two_unit_jan = parse(Float64, jan_row[1,10])
    threefour_unit_jan = parse(Float64, jan_row[1,13])
    total_jan = one_unit_jan + two_unit_jan + threefour_unit_jan

    one_unit_feb = parse(Float64, feb_row[1,7])
    two_unit_feb = parse(Float64, feb_row[1,10])
    threefour_unit_feb = parse(Float64, feb_row[1,13])
    total_feb = one_unit_feb + two_unit_feb + threefour_unit_feb
    missing_vec[i] = total_feb - total_jan
end

permits = fill(NaN, num_files, length(states))
permits[170, :] = missing_vec

for (t, file) in enumerate(readdir(save_dir))
    if file == "st1502c.txt" # 170
        continue
    end
    df = CSV.read(save_dir * file, DataFrame; header=1)
    for (i, state) in enumerate(states)
        state_row = df[(.!ismissing.(df.State)) .& (df.State .== state), :]
        one_unit = parse(Float64, state_row[1,7])
        two_unit = parse(Float64, state_row[1,10])
        threefour_unit = parse(Float64, state_row[1,13])
        permits[t, i] = one_unit + two_unit + threefour_unit
    end
end

save(datadir("updated_states/permits.jld2"), Dict("res" => permits))

# Some notes
# North Dakota in Jan/Feb of 2015 experienced actually negative housing permits
# Could be because of some corrections in the census data?

