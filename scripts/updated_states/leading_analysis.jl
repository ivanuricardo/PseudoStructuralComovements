using DrWatson
@quickactivate :PseudoStructuralComovements
using SparseArrays, MatrixAutoRegressions
Random.seed!(20260203)
include(projectdir("scripts/updated_states/helpers.jl"))

loaded_data = load(datadir("updated_states/leading_data.jld2"))

li_data = loaded_data["li_data"]
iowa_data = vcat(li_data[1:7:(end-1), :], li_data[end, :]')
illinois_data = vcat(li_data[2:7:(end-1), :], li_data[end, :]')
indiana_data = vcat(li_data[3:7:(end-1), :], li_data[end, :]')
minnesota_data = vcat(li_data[4:7:(end-1), :], li_data[end, :]')
northdakota_data = vcat(li_data[5:7:(end-1), :], li_data[end, :]')
ohio_data = vcat(li_data[6:7:(end-1), :], li_data[end, :]')
southdakota_data = vcat(li_data[7:7:(end-1), :], li_data[end, :]')

Plots.plot(demean_standardize(indiana_data[1:2,:]'))
Plots.plot!(demean_standardize(indiana_data[4,:]))

# Can not perform inference on the leading indicators as they take as input an
# estimated factor

iowa_var = restricted_VAR(iowa_data)
illinois_var = restricted_VAR(illinois_data)
indiana_var = restricted_VAR(indiana_data)
minnesota_var = restricted_VAR(minnesota_data)
northdakota_var = restricted_VAR(northdakota_data)
ohio_var = restricted_VAR(ohio_data)
southdakota_var = restricted_VAR(southdakota_data)


function forecast_CI(Y, A1, A2, A3, A4; steps=6)
    n, obs = size(Y)
    p = 4  # max lag in VAR

    Y_hist = Y[:, end-p+1:end]  # size: p × n

    Y_next = zeros(n)
    for h in 1:steps

        # CI equation: all 4 lags
        Y_next[1] = A1[1, :]' * Y_hist[:, end] +
                    A2[1, :]' * Y_hist[:, end-1] +
                    A3[1, :]' * Y_hist[:, end-2] +
                    A4[1, :]' * Y_hist[:, end-3]

        # X equations: only lag 1
        for i in 2:n
            Y_next[i] = A1[i, :]' * Y_hist[:, end]
        end

        # Append to history for next step
        Y_hist = hcat(Y_hist[:, 2:end], Y_next)
    end

    return Y_next[1]  # return CI forecast at horizon
end
