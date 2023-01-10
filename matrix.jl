function create_matrix_B(m::Int, h::Float64)::Matrix{Float64}

    B::Matrix{Float64} = zeros(m, m)

    f = x -> 1 / (h^2)

    B[1, 1] = -integral(f, 0, h * 2)
    B[1, 2] = integral(f, h, h * 2)
    B[m, m-1] = integral(f, h * m, h * (m + 1))
    B[m, m] = -integral(f, h * (m - 1), h * (m + 1))

    for i in 2:m-1
        B[i, i] = -integral(f, h * (i - 1), h * (i + 1))
        B[i, i+1] = integral(f, h * i, h * (i + 1))
        B[i, i-1] = integral(f, h * i, h * (i + 1))
        # println("1 = ", B[i, i], " 2 = ", B[i, i+1], " 0 = ", B[i, i-1])
    end

    return B
end

function create_matrix_L(G::Float64, a::Int, b::Int, e_fun, m::Int, h::Float64)

    L::Vector{Float64} = zeros(m)

    for j in 1:m
        if h * j > 1 && h * (j + 1) < 2
            L[j] = 4 * pi * G * 2 * integral(x -> e_fun(x, j, h), h * j, h * (j + 1), 2)
        end
        # println("hj = ", h * j, "     L[j] = ", L[j], "     hj+1 = ", h * (j + 1))
        # L[j] = 4 * pi * G * 2 * integral(h * j, h * (j + 1), 2)
    end

    return L
end