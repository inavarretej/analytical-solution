# Uncertainty quantification via finite-difference Hessian (Laplace approximation).
#
# Usage: evaluate the Hessian of χ² in log-parameter space at the optimum,
# then invert to get the covariance matrix for multiplicative 1σ uncertainties.

function hessian_fd(f, x; h=1e-3)
    n  = length(x)
    H  = zeros(n, n)
    f0 = f(x)
    for i in 1:n, j in i:n
        if i == j
            xp = copy(x); xp[i] += h
            xm = copy(x); xm[i] -= h
            H[i, i] = (f(xp) - 2*f0 + f(xm)) / h^2
        else
            xpp = copy(x); xpp[i] += h; xpp[j] += h
            xpm = copy(x); xpm[i] += h; xpm[j] -= h
            xmp = copy(x); xmp[i] -= h; xmp[j] += h
            xmm = copy(x); xmm[i] -= h; xmm[j] -= h
            H[i, j] = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4*h^2)
            H[j, i] = H[i, j]
        end
    end
    return H
end
