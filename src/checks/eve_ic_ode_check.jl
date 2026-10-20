using Printf
# One Fourier mode. Plate: t = K k^2 (u - g)  (basal traction balances plate stress, dike as eigen-displacement g)
# Maxwell half-space: (1/mu) dt/dt + t/eta = -2|k| du/dt    =>  du/dt = -( (K k^2/mu) dg/dt*(-1) ... ) solved below
K=40e9*10e3/(1-0.25^2); μ=16e9; η=5e18; k=2π/30e3; g0=1.0
ℓ=K/(2μ); a=K/(2η); τR=η/μ; α=a*k/(1+ℓ*k)
ε=1e-6*τR                                  # opening ramp duration
g(t)= t<=0 ? 0.0 : (t>=ε ? g0 : g0*t/ε); gd(t)= (0<t<ε) ? g0/ε : 0.0
# eliminate t: (K k^2/μ)(du - dg) + (K k^2/η)(u - g) = -2k du   (per unit time)
du(u,t)= ((K*k^2/μ)*gd(t) - (K*k^2/η)*(u-g(t))) / (K*k^2/μ + 2k)
u=0.0; t=0.0
function step!(u,t,h)  # RK4
  k1=du(u,t); k2=du(u+h/2*k1,t+h/2); k3=du(u+h/2*k2,t+h/2); k4=du(u+h*k3,t+h)
  u+h/6*(k1+2k2+2k3+k4)
end
n=2000; for i in 1:n; global u=step!(u,t,ε/n); global t+=ε/n; end
@printf("ℓk = %.3f\n", ℓ*k)
@printf("t=0+ : ODE u = %.6f   corrected ℓk/(1+ℓk) = %.6f   paper = 0\n", u, ℓ*k/(1+ℓ*k))
for tt in (0.5τR, 2τR, 10τR)
  h=(tt-t)/20000; for i in 1:20000; global u=step!(u,t,h); global t+=h; end
  @printf("t=%4.1f τR: ODE u = %.6f   corrected 1-e^{-αt}/(1+ℓk) = %.6f   paper 1-e^{-αt} = %.6f\n", tt/τR, u, 1-exp(-α*tt)/(1+ℓ*k), 1-exp(-α*tt))
end
