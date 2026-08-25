using FFTW

# Fourier convention matching the main text: fhat(k) = int f(x) e^{-ikx} dx.
# Hilbert transform symbol: F{Hf}(k) = -i*sign(k)*fhat(k).
# Derivative symbol:        F{f'}(k) = i*k*fhat(k).
# So F{d/dx[Hf]}(k) = (i*k)*(-i*sign(k))*fhat(k) = k*sign(k)*fhat(k) = |k|*fhat(k),
# which is exactly the symbol of (-Delta)^{1/2}. This script checks the claim
# on a discretized grid via FFT, for a smooth compactly-supported-in-practice
# test function (a Gaussian), and cross-checks against direct multiplication
# by |k| in Fourier space.

N  = 2^14
L  = 80.0
dx = L / N
x  = (-N÷2:N÷2-1) .* dx
k  = 2π .* fftfreq(N, 1/dx)

f = exp.(-x.^2)               # Gaussian test function
fhat = fft(f) .* dx            # approximate continuous FT via DFT

# Direct route: multiply by |k| in Fourier space, invert.
direct = real.(ifft(abs.(k) .* fhat) ./ dx)

# Hilbert-transform route: multiply by -i*sign(k), invert, differentiate via
# a second Fourier multiplication by i*k.
sgnk = sign.(k)
sgnk[k .== 0] .= 0.0
Hfhat = (-1im .* sgnk) .* fhat
dHfhat = (1im .* k) .* Hfhat
route2 = real.(ifft(dHfhat) ./ dx)

err = maximum(abs.(direct .- route2))
println("Max abs difference between |k|-route and d/dx[Hilbert]-route: ", err)
@assert err < 1e-6 "Hilbert transform / derivative symbol identity failed numerically"

println("All checks passed: (-Delta)^(1/2) = d/dx compose Hilbert transform, confirmed via FFT.")
