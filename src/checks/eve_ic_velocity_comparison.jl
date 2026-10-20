using CairoMakie, Printf, QuadGK
include(joinpath(pwd(),"src/models/viscous.jl")); include(joinpath(pwd(),"src/models/maxwell.jl"))
yr=3600.0*24*365.25; H=10e3; tR(a)=H/(a*(1-0.25))
# corrected post-diking RSS velocity: extra factor 1/(1+ℓk) (instantaneous elastic response removed from the relaxing amplitude)
function v_corr(x,τ,b,V,a,τR)
  T=b/V; ℓ=a*τR
  f(k)=(α=a*k/(1+ℓ*k); sin(k*x)/(1+ℓ*k)^2*exp(-α*τ)/-expm1(-α*T))
  (b*a/pi)*_segmented_quadgk(f,1e4/ℓ,x)
end
cases=[("Asal EVE fit",6.952e-5,0.642,13.0,19.0,12e3),("Krafla EVE fit",1.188e-4,1.599,18.0,14.5,150e3),("Fig. 8 setup, τ = 0.2T",nothing,nothing,nothing,nothing,nothing)]
fig=Figure(size=(1800,560),fontsize=20)
for (j,(name,a,b,Vmm,τy,xm)) in enumerate(cases)
  if a===nothing   # nondimensional-like illustration as in Fig. 8: T = 5 τR
    a=1e-4; τR=tR(a); T=5τR; b=1.0; V=b/T; τ=0.2T; xm=4.5*a*T; unit=T/b; ul="v / V"; xu=a*T; xl="x / aT"
  else
    V=Vmm/1000/yr; τ=τy*yr; τR=tR(a); unit=yr*1000; ul="v (mm/yr)"; xu=1e3; xl="x (km)"
  end
  xs=collect(range(-xm,xm,length=600))
  vp=[get_v_css_maxwell_split(x,τ,b,V,a,τR) for x in xs].*unit
  vc=[v_corr(x,τ,b,V,a,τR) for x in xs].*unit
  ve=[get_v_css(x,τ,b,V,a) for x in xs].*unit
  ax=Axis(fig[1,j],title=name,xlabel=xl,ylabel=ul)
  lines!(ax,xs./xu,ve,color=:purple,label="EV RSS (same a, b)")
  lines!(ax,xs./xu,vp,color=:navy,label="EVE, paper (no elastic step at t=0)")
  lines!(ax,xs./xu,vc,color=:darkorange,linewidth=3,label="EVE, corrected IC")
  j==1 && axislegend(ax,position=:rb,labelsize=15)
  @printf("%s: far field paper %.3f corrected %.3f EV %.3f\n",name,vp[end],vc[end],ve[end])
end
save(ARGS[1]*"/eve_ic_check.png",fig)
