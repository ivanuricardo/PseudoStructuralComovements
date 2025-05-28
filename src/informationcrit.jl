
function system_parameters(dimvals, ranks)
    first_term = ranks[1] * (dimvals[1] - ranks[1])
    second_term = ranks[2] * (dimvals[2] - ranks[2])
    third_term = ranks[1] * dimvals[1] + ranks[2] * dimvals[2]
    return first_term + second_term + third_term - 1
end

aic(ll::Real, numpars::Int) = -2 * ll + (2 * numpars)
bic(ll::Real, numpars::Int, obs::Int) = -2 * ll + (numpars * log(obs))
hqc(ll::Real, numpars::Int, obs::Int) = -2 * ll + (numpars * 2 * log(log(obs)))

function rank_selection(data, dimvals; iters=500)

    cen_data = data .- mean(data, dims=2)
    obs = size(cen_data, 2)
    ictable = fill(NaN, 5, prod(dimvals))
    rank_grid = collect(Iterators.product(1:dimvals[1], 1:dimvals[2]))

    for i = 1:prod(dimvals)
        #=Threads.@threads for i = 1:prod(dimvals)=#
        selected_rank = collect(rank_grid[i])
        num_parameters = system_parameters(dimvals, selected_rank)
        reg = comovement_reg(cen_data, dimvals, selected_rank; iters=iters)
        ll = -reg.res.minimum
        ictable[1, i] = aic(ll, num_parameters)
        ictable[2, i] = bic(ll, num_parameters, obs)
        ictable[3, i] = hqc(ll, num_parameters, obs)
        ictable[4, i] = selected_rank[1]
        ictable[5, i] = selected_rank[2]
    end

    aicvec = argmin(ictable[1, :])
    aic_sel = ictable[4:end, aicvec]
    bicvec = argmin(ictable[2, :])
    bic_sel = Int.(ictable[4:end, bicvec])
    hqcvec = argmin(ictable[3, :])
    hqc_sel = Int.(ictable[4:end, hqcvec])

    return (; aic_sel, bic_sel, hqc_sel, ictable)

end

function sim_stats(selectedvals::AbstractMatrix, correctval::AbstractVector, sims::Int)
    avgval = mean(selectedvals, dims=2)
    stdval = std(selectedvals, dims=2)
    numvals = size(selectedvals, 1)
    freqcorrect = fill(NaN, numvals)
    freqhigh = fill(NaN, numvals)
    freqlow = fill(NaN, numvals)

    for i = 1:numvals
        cval = correctval[i]
        freqcorrect[i] = count(x -> x == cval, selectedvals[i, :]) / sims
        freqhigh[i] = count(x -> x > cval, selectedvals[i, :]) / sims
        freqlow[i] = count(x -> x < cval, selectedvals[i, :]) / sims
    end

    return (; avgval, stdval, freqcorrect, freqhigh, freqlow)
end
