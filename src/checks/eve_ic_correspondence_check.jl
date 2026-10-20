# Textbook correspondence principle, one Fourier mode k, checked by numerical Laplace inversion.
# 1) Elastic problem: plate (K = EH/(1-ν²)) on an elastic half-space (shear modulus μ), dike opened by g0:
#      plate displacement d_el = g0 * ℓk/(1+ℓk),  ℓ = K/(2μ)
# 2) Correspondence: step load g0/s, and μ -> μ(s) = s·η̃(s) = s η/(1+τR s)
# 3) Invert numerically (Gaver–Stehfest) and compare with
#      corrected  d(t) = g0 [1 - e^{-αt}/(1+ℓ_M k)]   and   paper  d(t) = g0 [1 - e^{-αt}]
using Printf
setprecision(BigFloat, 200)
K = big(40e9) * 10e3 / (1 - 0.25^2); μm = big(16e9); η = big(5e18); k = 2π / big(30e3); g0 = big(1)
τR = η / μm; a = K / (2η); ℓM = K / (2μm); α = a * k / (1 + ℓM * k)
d_el(μ) = g0 * (K / (2μ)) * k / (1 + (K / (2μ)) * k)
d_laplace(s) = d_el(s * η / (1 + τR * s)) / s
function stehfest(F, t; N = 24)
    V(i) = sum(big(j)^(N ÷ 2) * factorial(big(2j)) /
               (factorial(big(N ÷ 2 - j)) * factorial(big(j)) * factorial(big(j - 1)) *
                factorial(big(i - j)) * factorial(big(2j - i)))
               for j in fld(i + 1, 2):min(i, N ÷ 2)) * (-1)^(N ÷ 2 + i)
    ln2t = log(big(2)) / t
    ln2t * sum(V(i) * F(i * ln2t) for i in 1:N)
end
@printf("ℓ_M k = %.3f\n", Float64(ℓM * k))
for f in (0.01, 0.5, 2.0, 10.0)
    t = f * τR
    @printf("t = %5.2f τR:  correspondence (Stehfest) = %.6f   corrected = %.6f   paper = %.6f\n", f,
            Float64(stehfest(d_laplace, t)), Float64(g0 * (1 - exp(-α * t) / (1 + ℓM * k))), Float64(g0 * (1 - exp(-α * t))))
end
