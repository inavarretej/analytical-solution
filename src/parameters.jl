# Physical constants for the Asal–Ghoubbet rifting analysis (Navarrete & Olive, 2026)

# Elastic moduli — Heap & Violay (2021, Bull. Volcanol. 83:33)
# Volcanic rock, porosity ϕ ≈ 0.05–0.10
lambda_L = 1.6e10      # Pa — Lamé first parameter
mu_L     = 1.6e10      # Pa — shear modulus
E_L      = mu_L * (3*lambda_L + 2*mu_L) / (lambda_L + mu_L)   # Young's modulus ≈ 40 GPa
nu_L     = lambda_L / (2*(lambda_L + mu_L))                     # Poisson ratio ≈ 0.25

# Asthenosphere — shared-modulus closure (paper Eq. mu_shared): the
# half-space is assumed to share the elastic plate's shear modulus,
# mu_A := mu_L, rather than an independently estimated mantle value.
mu_A  = mu_L           # Pa — mantle elastic shear modulus (= mu_L, shared closure)
# Viscosity — Nooner et al. (2009, GRL 36:L21313), Dabbahu segment, Afar:
# best-fit viscoelastic-relaxation viscosity for the East African Rift
# (same broader rift system as Natron), bracketing the range also seen at
# Krafla (Hofton & Foulger 1996) and Asal (this paper's own fit).
eta_a = 5.2e18         # Pa·s — asthenospheric viscosity

# Geometry
H = 10e3               # m — elastic plate thickness
h = 15e3               # m — diffusion depth scale (Foulger model)

# Time scales
yr    = 3600.0 * 24 * 365
tau   = 11 * yr         # s — observation epoch (11 yr post-event)
tau_8 = 8  * yr         # s — earlier epoch for differencing (8 yr post-event)
v     = 1 / (100 * yr)  # m/s — reference velocity

# Derived compound parameters
M     = 4 * mu_L * (lambda_L + mu_L) / (lambda_L + 2 * mu_L)
kappa = M * H * h / eta_a                         # m²/s — Foulger diffusivity
a_nom = (H * E_L) / (2 * eta_a * (1 - nu_L^2))   # m/s  — characteristic velocity
tau_R = eta_a / mu_A                               # s    — Maxwell relaxation time ≈ 10.3 yr
b_nom = 1.0                                        # m    — reference dike half-width
