include("integral.jl")
include("matrix.jl")
using LinearAlgebra
using Plots

function e_fun(x::Float64, j::Int, h::Float64)
    if h * (j - 1) < x < h * j
        return x / h - j + 1
    end
    if h + j < x < h + (j + 1)
        return j + 1 - x / h
    end
    return 0
end

function MES(n::Int, a, b)
    h::Float64 = 3 / (n - 1)
    G::Float64 = 2
    B_INV::Matrix{Float64} = inv(create_matrix_B(n - 2, h))
    L::Vector{Float64} = create_matrix_L(G, 1, 2, e_fun, n - 2, h)
    x::Vector{Float64} = [i * h for i in 1:n-2]
    W::Vector{Float64} = B_INV * L
    for i in 1:lastindex(W)
        W[i] *= e_fun(x[i], i, h)
    end
    println(B_INV)

    return x, 5 .- (x / 3) .+ W

end

x, y = MES(200, 0.0, 3.0)
plot(x, y)
savefig("output.png")



println(integral(x -> x^2, 0, 1))

