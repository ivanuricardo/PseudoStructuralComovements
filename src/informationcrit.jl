
function system_parameters(dimvals, ranks; p=1)
    first_term = ranks[1] * (dimvals[1] - ranks[1])
    second_term = ranks[2] * (dimvals[2] - ranks[2])
    third_term = p * ranks[1] * dimvals[1] + p * ranks[2] * dimvals[2] - p
    num_ll1 = Int(dimvals[1] * (dimvals[1] + 1) / 2) - 1
    num_ll2 = Int(dimvals[2] * (dimvals[2] + 1) / 2)

    return first_term + second_term + num_ll1 + num_ll2 + third_term
end

aic(ll::Real, numpars::Int) = -2 * ll + (2 * numpars)
bic(ll::Real, numpars::Int, obs::Int) = -2 * ll + (numpars * log(obs))
function ebic(ll::Real, dimvals::AbstractVector, ranks::AbstractVector, obs::Int)
    n1, n2 = dimvals
    r1, r2 = ranks
    firstdim = log(obs * n2) * r1 * (2 * n1 - r1)
    seconddim = log(obs * n1) * r2 * (2 * n2 - r2)
    return -2 * ll + firstdim + seconddim
end
hqc(ll::Real, numpars::Int, obs::Int) = -2 * ll + (numpars * 2 * log(log(obs)))

function rank_selection(data, dimvals; iters=1000, pmax=1)

    obs = size(data, 2)
    ictable = fill(NaN, 6, prod(dimvals) * pmax)
    rank_grid = collect(Iterators.product(1:dimvals[1], 1:dimvals[2], 1:pmax))

    # for i = 1:(prod(dimvals) * pmax)
    @showprogress Threads.@threads for i = 1:(prod(dimvals) * pmax)
        selected_rank_lags = collect(rank_grid[i])
        selected_rank = selected_rank_lags[1:2]
        selected_lag = selected_rank_lags[3]

        num_parameters = system_parameters(dimvals, selected_rank; p=selected_lag)
        # if selected_lag == 1
        #     reg = rrmar(data, dimvals, selected_rank)
        #     ll = -reg.ll
        # else
        reg = comovement_reg(data, dimvals, selected_rank; iters=iters, p=selected_lag)
        ll = -reg.res.minimum
        # end
        ictable[1, i] = aic(ll, num_parameters)
        ictable[2, i] = bic(ll, num_parameters, obs)
        ictable[3, i] = ebic(ll, dimvals, selected_rank, obs)
        ictable[4, i] = selected_rank[1]
        ictable[5, i] = selected_rank[2]
        ictable[6, i] = selected_lag
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
    mad = mean(abs.(selectedvals .- correctval), dims=2)
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

    return (; mad, stdval, freqcorrect, freqhigh, freqlow)
end
