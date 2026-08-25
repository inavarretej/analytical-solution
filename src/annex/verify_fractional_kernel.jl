using SpecialFunctions, QuadGK

# ── (a) C(N,s) at N=1, s=1/2 ────────────────────────────────────────────────
C(N, s) = s * 4.0^s / π^(N/2) * gamma(s + N/2) / gamma(1 - s)

C1half = C(1, 0.5)
println("C(1, 1/2)  = ", C1half)
println("1/pi       = ", 1/π)
println("difference = ", abs(C1half - 1/π))
@assert abs(C1half - 1/π) < 1e-12 "C(1,1/2) does not match 1/pi"

# ── (b) P.V. singular-difference integral vs. Hadamard finite part ─────────
# Test function: Gaussian bump, smooth and rapidly decaying.
f(y) = exp(-y^2)
x0 = 0.7  # evaluation point, away from the bump center to avoid symmetry
           # accidentally hiding a sign error

# P.V. of [f(x0)-f(y)]/(x0-y)^2 over y, computed as a genuine principal value
# by symmetric exclusion of a shrinking interval around y = x0.
function pv_diff(f, x0; eps=1e-6, lim=40.0)
    g(y) = (f(x0) - f(y)) / (x0 - y)^2
    left,  _ = quadgk(g, -lim, x0 - eps, rtol=1e-10)
    right, _ = quadgk(g, x0 + eps, lim, rtol=1e-10)
    return left + right
end

# Hadamard finite part of f(y)/(x0-y)^2 over y: subtract the divergent part
# before taking eps -> 0. The analytic infinite-domain divergence is
# 2*f(x0)/eps, but this quadrature is truncated at +-lim, so the raw
# integral only contains int_{eps<=|z|<=lim} dz/z^2 = 2/eps - 2/lim of that
# term, not the full 2/eps. Subtracting the untruncated 2/eps leaves a
# systematic O(f(x0)/lim) residual (e.g. ~0.03 at lim=40) that has nothing
# to do with eps -> 0 convergence and does not shrink as eps shrinks -- it
# is a finite-domain bookkeeping mismatch. Subtract the truncation-matched
# constant instead so the two computations are regularizing over the same
# domain.
function hadamard_fp(f, x0; eps=1e-6, lim=40.0)
    g(y) = f(y) / (x0 - y)^2
    left,  _ = quadgk(g, -lim, x0 - eps, rtol=1e-10)
    right, _ = quadgk(g, x0 + eps, lim, rtol=1e-10)
    return (left + right) - f(x0)*(2/eps - 2/lim)
end

pv_val  = pv_diff(f, x0; eps=1e-5)
hfp_val = hadamard_fp(f, x0; eps=1e-5)

println("\nP.V. integral      = ", pv_val)
println("-Hfp integral      = ", -hfp_val)
println("difference         = ", abs(pv_val - (-hfp_val)))
@assert abs(pv_val - (-hfp_val)) < 1e-3 "P.V. and -Hfp do not agree"

println("\nAll checks passed.")
