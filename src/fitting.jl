# Chi-squared misfit functions, parallel evaluator, and Latin hypercube sampler.
# All chi² functions take observed data explicitly — wrap with closures in the
# calling script/notebook to bind data once and pass single-argument callables
# to pmap_chi2 / optimize.
#
# Dependencies (must be loaded before include):
#   QuasiMonteCarlo (sample, LatinHypercubeSample)
#   models.jl (diff_step, diff_foulger, diff_viscoelastic)

# ── Parallel evaluator ─────────────────────────────────────────────────────────

function pmap_chi2(f, params)
    n       = length(params)
    results = Vector{Float64}(undef, n)
    Threads.@threads for i in 1:n
        results[i] = f(params[i])
    end
    return results
end

# ── Latin hypercube sampler (log-uniform) ──────────────────────────────────────

function lhs_log(lb, ub, n)
    plan = QuasiMonteCarlo.sample(n, log10.(lb), log10.(ub), LatinHypercubeSample())
    return [10 .^ plan[:, i] for i in 1:n]
end

# ── 1-parameter misfits (b fixed) ─────────────────────────────────────────────

function chi2_step(a_test, obs_x, obs_u, obs_sig, tau_i, tau_f, b)
    u_pred = diff_step(obs_x, tau_i, tau_f, b, a_test)
    return sum(((obs_u .- u_pred) ./ obs_sig).^2)
end

function chi2_foulger(k_test, obs_x, obs_u, obs_sig, tau_i, tau_f, b)
    u_pred = diff_foulger(obs_x, tau_i, tau_f, b, k_test)
    return sum(((obs_u .- u_pred) ./ obs_sig).^2)
end

function chi2_viscoelastic_1p(a_test, obs_x, obs_u, obs_sig, tau_i, tau_f, b, tau_R)
    u_pred = diff_viscoelastic(obs_x, tau_i, tau_f, b, a_test, tau_R)
    return sum(((obs_u .- u_pred) ./ obs_sig).^2)
end

# ── 2-parameter misfits (a or κ, and b jointly) ───────────────────────────────

function chi2_step_2p(params, obs_x, obs_u, obs_sig, tau_i, tau_f)
    a_val, b_val = params
    (a_val <= 0 || b_val <= 0) && return 1e10
    return sum(((obs_u .- diff_step(obs_x, tau_i, tau_f, b_val, a_val)) ./ obs_sig).^2)
end

function chi2_foulger_2p(params, obs_x, obs_u, obs_sig, tau_i, tau_f)
    k_val, b_val = params
    (k_val <= 0 || b_val <= 0) && return 1e10
    return sum(((obs_u .- diff_foulger(obs_x, tau_i, tau_f, b_val, k_val)) ./ obs_sig).^2)
end

function chi2_viscoelastic_2p(params, obs_x, obs_u, obs_sig, tau_i, tau_f, tau_R)
    a_val, b_val = params
    (a_val <= 0 || b_val <= 0) && return 1e10
    return sum(((obs_u .- diff_viscoelastic(obs_x, tau_i, tau_f, b_val, a_val, tau_R)) ./ obs_sig).^2)
end

# ── 3-parameter Maxwell misfit (a, b, τ_R all free) ───────────────────────────

function chi2_viscoelastic_3p(params, obs_x, obs_u, obs_sig, tau_i, tau_f)
    a_val, b_val, tauR_val = params
    (a_val <= 0 || b_val <= 0 || tauR_val <= 0) && return 1e10
    u_pred = diff_viscoelastic(obs_x, tau_i, tau_f, b_val, a_val, tauR_val)
    return sum(((obs_u .- u_pred) ./ obs_sig).^2)
end
