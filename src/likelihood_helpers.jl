
function blockvec(A::AbstractMatrix, dimvals, ranks)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    A1 = A[1:N1_r1, 1:N2_r2]
    A2 = A[(N1_r1+1):end, 1:N2_r2]
    A3 = A[1:N1_r1, (N2_r2+1):end]
    A4 = A[(N1_r1+1):end, (N2_r2+1):end]
    return vcat(vec(A1), vec(A2), vec(A3), vec(A4))
end

# Creates the permutation matrix P such that
# Ω P' P y_t works out
function perm_matrix(dimvals, ranks)
    seq = reshape(1:prod(dimvals), dimvals[1], dimvals[2])
    res = I(prod(dimvals))
    return res[:, blockvec(seq, dimvals, ranks)]'
end

top_left(N1_r1, N2_r2) = I(N1_r1 * N2_r2)

function right1(delta_star, N2_r2)
    gamma_eye = I(N2_r2)
    tr_kron = kron(gamma_eye, delta_star')
    return tr_kron
end

function right2(gamma_star, N1_r1)
    delta_eye = I(N1_r1)
    tr_kron = kron(gamma_star', delta_eye)
    return tr_kron
end

function right3(delta_star, gamma_star)
    tr_kron = kron(gamma_star', delta_star')
    return tr_kron
end

function create_second_row(delta_star, ranks, N1_r1, N2_r2)
    zero_dim = N1_r1 * N2_r2 + ranks[1] * N2_r2
    left_z = zeros(ranks[2] * N1_r1, zero_dim)
    right_eye = I(ranks[2] * N1_r1)
    right_kron = kron(I(ranks[2]), delta_star')
    return [left_z right_eye right_kron]
end
function create_third_row(gamma_star, ranks, N1_r1, N2_r2)
    left_z = zeros(ranks[1] * N2_r2, N1_r1 * N2_r2)
    left_eye = I(ranks[1] * (N2_r2))
    right_z = zeros(ranks[1] * N2_r2, ranks[2] * N1_r1)
    right_kron = kron(gamma_star', I(ranks[1]))
    return [left_z left_eye right_z right_kron]
end

function create_fourth_row(dimvals, ranks)
    r1_r2 = ranks[1] * ranks[2]
    left_zeros = zeros(r1_r2, prod(dimvals) - r1_r2)
    right_eye = I(r1_r2)
    return [left_zeros right_eye]
end

function create_omega(delta_star, gamma_star, dimvals, ranks)
    N1_r1 = dimvals[1] - ranks[1]
    N2_r2 = dimvals[2] - ranks[2]

    t1 = top_left(N1_r1, N2_r2)
    t2 = right1(delta_star, N2_r2)
    t3 = right2(gamma_star, N1_r1)
    t4 = right3(delta_star, gamma_star)

    first_row = hcat(t1, t2, t3, t4)
    second_row = create_second_row(delta_star, ranks, N1_r1, N2_r2)
    third_row = create_third_row(gamma_star, ranks, N1_r1, N2_r2)
    fourth_row = create_fourth_row(dimvals, ranks)
    omega = vcat(first_row, second_row, third_row, fourth_row)

    return omega
end

function create_pi(u3, u4, dimvals, ranks; p=1)
    pdims = prod(dimvals)
    r = prod(ranks)
    if p == 1
        top_pi = zeros(pdims - r, pdims)
        bottom_pi = kron(u4, u3)'
        return vcat(top_pi, bottom_pi)
    end

    m = pdims - r
    t_u3 = eltype(u3)
    pi_result = zeros(t_u3, pdims, pdims * p)
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

