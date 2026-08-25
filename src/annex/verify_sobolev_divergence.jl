using QuadGK

# u(x) = H(-x), the Heaviside step (the model's dike-opening initial
# condition up to the scale factor b: u(x,0+) = b*H(-x)). The Gagliardo
# seminorm for H^{1/2}(R) is
#   [u]^2 = int int |u(x)-u(y)|^2 / (x-y)^2  dx dy.
# Restricting to opposite-sign pairs (the only pairs where u(x) != u(y)):
#   [u]^2 = 2 * I,   I = int_0^inf int_0^inf dt dx / (x+t)^2
# (the prefactor is the shared divergent core; it does not depend on
# whether the step has size 1 (Heaviside) or 2 (sign) -- only an overall
# constant multiplies it, and that constant does not affect whether I
# diverges). The inner integral over t is exact: int_0^inf dt/(x+t)^2 = 1/x.
# So I = int_0^inf dx / x, which diverges logarithmically at x -> 0.
# We confirm this numerically by truncating the divergent end at delta and
# showing I(delta) grows like -log(delta) as delta -> 0.

inner(x) = quadgk(t -> 1/(x + t)^2, 0.0, 40.0, rtol=1e-10)[1]  # should be ~ 1/x
xs = [0.5, 0.05, 0.005]
println("Check inner integral ~ 1/x:")
for x in xs
    println("  x=$x  inner=$(inner(x))  1/x=$(1/x)")
end

function I_truncated(delta; upper=40.0)
    quadgk(x -> inner(x), delta, upper, rtol=1e-8)[1]
end

deltas = [1e-2, 1e-3, 1e-4, 1e-5]
vals = [I_truncated(d) for d in deltas]

println("\nTruncated I(delta):")
for (d, v) in zip(deltas, vals)
    println("  delta=$d  I=$v  -log(delta)=$(-log(d))")
end

# Confirm I(delta) - (-log(delta)) approaches a constant (i.e. I(delta) grows
# exactly like -log(delta) + O(1), not faster or slower).
offsets = vals .- (-log.(deltas))
println("\nI(delta) + log(delta)  (should converge to a constant):")
for (d, o) in zip(deltas, offsets)
    println("  delta=$d  offset=$o")
end
spread = maximum(offsets) - minimum(offsets)
println("\nSpread across four orders of magnitude in delta: ", spread)
@assert spread < 0.05 "I(delta) does not track -log(delta) as expected"

println("\nAll checks passed: the Gagliardo seminorm of H(-x) diverges logarithmically.")
