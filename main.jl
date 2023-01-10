include("integral.jl")
include("matrix.jl")
using LinearAlgebra
using Plots

function e_fun(x::Float64, j::Int, h::Float64)
    if h * (j - 1) < x < h * j
        return (x / h) - j + 1
    end
    if h * j < x < h * (j + 1)
        return j + 1 - (x / h)
    end
    return 0
end

function MES(n::Int)
    h::Float64 = 3 / (n - 1)
    G::Float64 = 2
    B = create_matrix_B(n - 2, h)
    B_INV::Matrix{Float64} = inv(B)
    L::Vector{Float64} = create_matrix_L(G, 1, 2, e_fun, n - 2, h)
    A::Vector{Float64} = B_INV * L

    w(x) = sum([e_fun(x, i, h) * A[i] for i in 1:(n-2)])

    return w
end

w = @time MES(100)

a::Float64, b::Float64 = 0.0, 3.0
p::Int = 800
dx::Float64 = (b - a) / p
x::Vector{Float64} = a:dx:b
y = w.(x) .+ 5 .- (x / 3)
plot(x, y)
savefig("output.png")


