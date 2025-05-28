
function both_blockvec(A::AbstractMatrix, dimvals, ranks)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    A1 = A[1:N1_r1, 1:N2_r2]
    A2 = A[(N1_r1+1):end, 1:N2_r2]
    A3 = A[1:N1_r1, (N2_r2+1):end]
    A4 = A[(N1_r1+1):end, (N2_r2+1):end]
    return vcat(vec(A1), vec(A2), vec(A3), vec(A4))
end

function both_perm_mat(dimvals, ranks)
    seq = reshape(1:prod(dimvals), dimvals[1], dimvals[2])
    res = I(prod(dimvals))
    return res[:, both_blockvec(seq, dimvals, ranks)]'
end

b_top_left(N1_r1, N2_r2) = I(N1_r1 * N2_r2)

function b_right1(delta_star, N2_r2)
    gamma_eye = I(N2_r2)
    tr_kron = kron(gamma_eye, delta_star')
    return tr_kron
end

function b_right2(gamma_star, N1_r1)
    delta_eye = I(N1_r1)
    tr_kron = kron(gamma_star', delta_eye)
    return tr_kron
end

function b_right3(delta_star, gamma_star)
    tr_kron = kron(gamma_star', delta_star')
    return tr_kron
end

function b_secondrow(delta_star, ranks, N1_r1, N2_r2)
    zero_dim = N1_r1 * N2_r2 + ranks[1] * N2_r2
    left_z = zeros(ranks[2] * N1_r1, zero_dim)
    right_eye = I(ranks[2] * N1_r1)
    right_kron = kron(I(ranks[2]), delta_star')
    return [left_z right_eye right_kron]
end
function b_thirdrow(gamma_star, ranks, N1_r1, N2_r2)
    left_z = zeros(ranks[1] * N2_r2, N1_r1 * N2_r2)
    left_eye = I(ranks[1] * (N2_r2))
    right_z = zeros(ranks[1] * N2_r2, ranks[2] * N1_r1)
    right_kron = kron(gamma_star', I(ranks[1]))
    return [left_z left_eye right_z right_kron]
end

function b_fourthrow(dimvals, ranks)
    r1_r2 = ranks[1] * ranks[2]
    left_zeros = zeros(r1_r2, prod(dimvals) - r1_r2)
    right_eye = I(r1_r2)
    return [left_zeros right_eye]
end

function omega_from_both(delta_star, gamma_star, dimvals, ranks)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    t1 = b_top_left(N1_r1, N2_r2)
    t2 = b_right1(delta_star, N2_r2)
    t3 = b_right2(gamma_star, N1_r1)
    t4 = b_right3(delta_star, gamma_star)

    first_row = hcat(t1, t2, t3, t4)
    second_row = b_secondrow(delta_star, ranks, N1_r1, N2_r2)
    third_row = b_thirdrow(gamma_star, ranks, N1_r1, N2_r2)
    fourth_row = b_fourthrow(dimvals, ranks)
    omega = vcat(first_row, second_row, third_row, fourth_row)

    return omega
end

function pi_from_both(u3, u4, dimvals, ranks; p=1)
    pdims = prod(dimvals)
    r = prod(ranks)
    if p == 1
        top_pi = zeros(pdims - r, pdims)
        bottom_pi = kron(u4, u3)'
        return vcat(top_pi, bottom_pi)
    end

    m = pdims - r
    pi_result = zeros(pdims, pdims * p)
    for k in 0:(p-1)
        col1 = k * pdims + 1
        col2 = (k + 1) * pdims

        i1, i2 = k * dimvals[1] + 1, (k + 1) * dimvals[1]
        j1, j2 = k * dimvals[2] + 1, (k + 1) * dimvals[2]

        bottom = kron(
            @view(u4[j1:j2, :]),
            @view(u3[i1:i2, :])
        )'

        pi_result[m+1:end, col1:col2] = bottom
    end
    return pi_result
end

function both_loglike(delta, gamma, u3, u4, sigma, resp, pred, dimvals, ranks; p=1)
    obs = size(resp, 2)

    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    rotate_u!(delta)
    delta_star = delta[(N1_r1+1):end, :]
    rotate_u!(gamma)
    gamma_star = gamma[(N2_r2+1):end, :]

    omega = omega_from_both(delta_star, gamma_star, dimvals, ranks)
    pi_mat = pi_from_both(u3, u4, dimvals, ranks)
    sse = 0.0

    if p == 1
        logdet_term = logdet(omega * sigma * omega')
        precision_matrix = inv(omega * sigma * omega')

        for i = 2:obs
            yt = resp[:, i]
            yt_m1 = pred[:, i]
            resid = omega * yt - pi_mat * yt_m1
            sse += dot(resid, precision_matrix * resid)
        end

        return 0.5 * (obs * logdet_term + sse)
    end

    omega_tilde, pi_tilde = make_companion(omega, pi_mat)
    logdet_term = logdet(omega_tilde * sigma * omega_tilde')
    precision_matrix = inv(omega_tilde * sigma * omega')
    for i = 2:obs
        yt = resp[:, i]
        yt_m1 = pred[:, i]
        resid = omega * yt - pi_mat * yt_m1
        sse += dot(resid, precision_matrix * resid)
    end




end
