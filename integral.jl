using FastGaussQuadrature

function integral(f, a, b, n=2)
    "Compute the integral using the quadrature Gauss"

    # change a, b to -1, 1
    scaling = (b - a) / 2
    shift = (b + a) / 2

    # compute the nodes and weights
    x, w = gausslegendre( n )
    return w' * f.(x .* scaling .+ shift)
end


