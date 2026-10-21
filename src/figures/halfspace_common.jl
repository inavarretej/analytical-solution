using CairoMakie, QuadGK, SpecialFunctions

include("../models/viscous.jl")
include("../models/maxwell.jl")
include("../parameters.jl")
include("../halfspace.jl")

mkpath("fig")

# a_nom comes from src/parameters.jl: a = H*E_L/[2*eta_a*(1-nu_L^2)], derived
# from the literature-sourced plate modulus (Heap & Violay 2021) and
# asthenospheric viscosity (Nooner et al. 2009, Afar), NOT an independent
# best-fit value -- so this figure is self-consistent with the EVE companion
# figure's tau_R = eta_a/mu_A under the shared-modulus closure.
b_nom = 0.662       # m

# Snapshot phases as multiples of tau_R (src/parameters.jl), not fixed yr
# values: since tau_EV and tau_R now coincide under the shared-modulus
# closure (mu_A := mu_L), expressing both this figure's and the EVE
# companion figure's snapshots as the same multiples of tau_R puts them on
# a directly comparable time axis.
tau_R_mult = [0.001, 0.25, 0.5, 1.0, 2.0, 4.9]
taus_yr = tau_R_mult .* (tau_R / yr)
taus    = taus_yr .* yr
labels  = ["$(m)·τ_R" for m in tau_R_mult]

x_span_m = 300e3    # +/- 300 km
z_top_m  = 20e3     # 20 km: shallow margin below the surface (also sets the room available for arrow length)
z_bot_m  = 300e3    # 300 km depth
x_prof = collect(range(-x_span_m, x_span_m, length = 400))

# Interior quiver keeps its own, narrower x-range (independent of x_span_m
# above): widening it in step with the surface-profile row would have halved
# the on-screen arrow length (DataAspect ties x and z to the same physical
# scale, and z_top_m/nz -- unrelated to x_span_m -- already caps arrow length
# to avoid overlapping the row above/below it).
x_quiver_span_m = 300e3

nx, nz = 15, 9
x_grid = collect(range(-x_quiver_span_m, x_quiver_span_m, length = nx))
z_grid = collect(range(z_top_m, z_bot_m, length = nz))
xs = Float64[]; zs = Float64[]
for zi in z_grid, xi in x_grid
    push!(xs, xi); push!(zs, zi)
end

# Exclude the grid point sitting right above the dike (x=0, shallowest z):
# the near-source field there is huge compared to the rest of the domain and
# was dominating both the color scale and the arrow-length scale.
origin_mask = isapprox.(xs, 0.0; atol = 1.0) .& isapprox.(zs, z_top_m; atol = 1.0)
keep  = .!origin_mask

# The far corners (max |x|, max z) excluded from the *color-range*
# computation only (compute_speed_min/compute_speed_max below), not from the
# plotted field itself: those points sit much farther from the dike than
# anything else sampled, so at the earliest phase they can pull a dynamic
# color_min down by orders of magnitude, stretching the whole log-scale
# range without carrying real information (the rest of the domain then
# bunches up near the top of that range).
corner_mask     = isapprox.(abs.(xs), x_quiver_span_m; atol = 1.0) .& isapprox.(zs, z_bot_m; atol = 1.0)
keep_for_range  = keep .& .!corner_mask

# Dense structured grid for the :contour interior panel (plot_model_rows!)
# -- much finer than the sparse quiver grid above, and used unmasked (no
# origin exclusion): the fixed color_min/color_max range already clamps the
# near-dike spike to the top color, same as out-of-range arrow colors do.
# Starts at z=0 (the surface itself), unlike the sparse quiver grid's
# z_top_m margin: at z=0, halfspace_vx/vz reduce to the surface profile
# V(x)/0 exactly (no singularity), and nx_dense even means x=0 is never a
# sampled column, so the one truly singular point (x=0, z=0) is never hit.
nx_dense, nz_dense = 400, 240
x_grid_dense = collect(range(-x_quiver_span_m, x_quiver_span_m, length = nx_dense))
z_grid_dense = collect(range(0.0, z_bot_m, length = nz_dense))

# Shared fixed color range for all four halfspace_{step,css}[_contour].jl
# figures: covers both models' actual data (STEP ~3e-4 to 0.08, CSS ~0.17 to
# 0.58, see keep_for_range above) with round bounds, on a log10 colorscale --
# used by every script instead of an autoscaled range so the color of a
# given |v| means the same thing across STEP/CSS and quiver/contour.
SHARED_COLOR_MIN = 1e-4
SHARED_COLOR_MAX = 1.5

n_rows = length(taus)  # number of phase snapshots; used by each script to size its Figure

# Filename tag for a generate_*_figure(; orientation, ...) call: empty for
# the default :vertical (keeps existing filenames, e.g. fig_halfspace_step.png,
# so the paper's \includegraphics references don't need updating), "_horizontal"
# otherwise.
orientation_suffix(orientation::Symbol) = orientation == :vertical ? "" : "_$(orientation)"

function compute_speed_max(taus, make_V, nondim)
    speed_max = 0.0
    for tau_i in taus
        V = make_V(tau_i)
        us = halfspace_vx.(xs, zs, Ref(V)) .* nondim
        ws = halfspace_vz.(xs, zs, Ref(V)) .* nondim
        speed_max = max(speed_max, maximum(sqrt.(us .^ 2 .+ ws .^ 2)[keep_for_range]))
    end
    return speed_max
end

function compute_speed_min(taus, make_V, nondim)
    speed_min = Inf
    for tau_i in taus
        V = make_V(tau_i)
        us = halfspace_vx.(xs, zs, Ref(V)) .* nondim
        ws = halfspace_vz.(xs, zs, Ref(V)) .* nondim
        speed_min = min(speed_min, minimum(sqrt.(us .^ 2 .+ ws .^ 2)[keep_for_range]))
    end
    return speed_min
end

# Dense-grid counterparts of compute_speed_max/compute_speed_min, for
# interior_panel = :contour: min/max over the SAME dense grid
# plot_model_rows! actually contours, rather than the sparse quiver grid --
# "the actual velocities" being plotted -- EXCLUDING a near-source box
# (|x|<exclude_extent[1], z<exclude_extent[2]; default matches the
# interior_panel zoom_extent default) from the min/max search only, not
# from the plotted field. Without this, the near-dike spike sets color_max
# for every phase from a single one of them (typically the earliest), so a
# log-scale range wide enough to include it leaves the other phases' own
# (much smaller) dynamic range flat/indistinguishable near the bottom.
function compute_dense_speed_max(taus, make_V, nondim; exclude_extent = (10e3, 10e3))
    x_max_excl, z_max_excl = exclude_extent
    speed_max = 0.0
    for tau_i in taus
        V = make_V(tau_i)
        for zi in z_grid_dense, xi in x_grid_dense
            (abs(xi) < x_max_excl && zi < z_max_excl) && continue
            us = halfspace_vx(xi, zi, V) * nondim
            ws = halfspace_vz(xi, zi, V) * nondim
            speed_max = max(speed_max, sqrt(us^2 + ws^2))
        end
    end
    return speed_max
end

function compute_dense_speed_min(taus, make_V, nondim; exclude_extent = (10e3, 10e3))
    x_max_excl, z_max_excl = exclude_extent
    speed_min = Inf
    for tau_i in taus
        V = make_V(tau_i)
        for zi in z_grid_dense, xi in x_grid_dense
            (abs(xi) < x_max_excl && zi < z_max_excl) && continue
            us = halfspace_vx(xi, zi, V) * nondim
            ws = halfspace_vz(xi, zi, V) * nondim
            speed_min = min(speed_min, sqrt(us^2 + ws^2))
        end
    end
    return speed_min
end

# Progress logging for long runs: timestamped, flushed immediately so a log
# file (stdout redirected) shows progress in real time.
using Dates, Serialization
fmt_hms(s) = (s = round(Int, s); string(s ÷ 3600, "h", lpad((s % 3600) ÷ 60, 2, '0'), "m", lpad(s % 60, 2, '0'), "s"))
function log_progress(msg)
    println(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "  ", msg)
    flush(stdout)
end

# Checkpointing of expensive grids (see cache_prefix in plot_model_rows!).
function cached_grid(f, prefix, name)
    prefix === nothing && return f()
    path = joinpath("cache", "$(prefix)_$(name).jls")
    if isfile(path)
        log_progress("$name loaded from $path")
        return deserialize(path)
    end
    g = f()
    mkpath("cache")
    serialize(path, g)
    log_progress("$name saved to $path")
    return g
end

# One model (STEP or CSS) per figure. `orientation = :vertical` lays phases
# out as rows (column 1 = surface velocity/displacement profile, column 2 =
# interior velocity quiver, one shared colorbar column); `:horizontal` swaps
# rows and columns, so phases run left-to-right and the profile/quiver panels
# each get their own row, with the shared colorbar in a trailing column.
#
# Each model has its own physically appropriate normalization, not a shared
# one: STEP distances by the plate thickness H, displacement by b, velocity
# by b/tau_R (tau_R being the Maxwell time set by the half-space viscosity
# and the elastic layer's modulus); CSS distance AND displacement by the
# fixed length aT (the RSS's own characteristic scale, matching
# Figure~fig:css_behavior's convention), velocity by the plate rate V=b/T,
# with phases expressed as fractions of the recurrence interval T rather
# than of tau_R.
function plot_model_rows!(fig, taus, labels, make_V, make_U;
        orientation::Symbol = :vertical, fontsize::Real = 45.8,
        interior_panel::Symbol = :quiver, zoom_extent = nothing, title, colorbar_label,
        vel_nondim, vel_symbol, disp_nondim, disp_symbol, u_lim, vel_lim,
        length_scale, length_symbol, guide_formula, color_min, color_max,
        ref_speed, colorscale = identity,
        hvx::Function = halfspace_vx, hvz::Function = halfspace_vz,
        guide_length_fn::Function = tau -> a_nom * tau,
        time_symbol::AbstractString = "τ",
        cache_prefix = nothing)
    # hvx/hvz: interior-field kernels used by the :contour panels (the EVE
    # scripts pass maxevals-capped versions; the uncapped defaults hang on the
    # nested Maxwell quadrature). guide_length_fn: position of the dashed
    # vertical guide in the profile panels (default L = a*tau).
    # cache_prefix: if given, each finished contour grid is saved to
    # cache/<prefix>_<name>.jls and reloaded on the next run (resume after an
    # interruption). Delete cache/<prefix>_* whenever the model or grid changes.

    # |v| on a structured grid, threaded over points (deterministic: same
    # values as a serial loop, only faster). Logs progress every ~2% with
    # elapsed time and ETA.
    function speed_grid(xg, zg, V; tag = "")
        sp = Matrix{Float64}(undef, length(xg), length(zg))
        n_tot = length(sp)
        n_done = Threads.Atomic{Int}(0)
        step = max(1, n_tot ÷ 50)
        t0 = time()
        Threads.@threads for idx in eachindex(sp)
            ix, iz = Tuple(CartesianIndices(sp)[idx])
            ux = hvx(xg[ix], zg[iz], V) * vel_nondim
            uz = hvz(xg[ix], zg[iz], V) * vel_nondim
            sp[ix, iz] = sqrt(ux^2 + uz^2)
            c = Threads.atomic_add!(n_done, 1) + 1
            if c % step == 0
                el = time() - t0
                log_progress("$tag  $(round(100 * c / n_tot, digits = 1))%  elapsed $(fmt_hms(el))  ETA $(fmt_hms(el * (n_tot - c) / c))")
            end
        end
        return sp
    end

    grid_or_cache(f, name) = cached_grid(f, cache_prefix, name)

    orientation in (:vertical, :horizontal) ||
        error("orientation must be :vertical or :horizontal, got $(repr(orientation))")
    interior_panel in (:quiver, :contour) ||
        error("interior_panel must be :quiver or :contour, got $(repr(interior_panel))")
    zoom_extent === nothing || interior_panel == :contour ||
        error("zoom_extent is only supported with interior_panel = :contour")

    # Title/legend/reference-text sizes scale with `fontsize` (the Figure's
    # own base size), keeping the ratios the layout was tuned at (base 45.8).
    legend_fontsize = fontsize * (41.63 / 45.8)
    title_fontsize  = fontsize * (54.12 / 45.8)
    ref_fontsize    = fontsize * (33.31 / 45.8)

    n = length(taus)
    # Number of panel columns (rows, under :horizontal) per phase: profile +
    # interior, plus a third zoom panel when zoom_extent is given.
    n_panels = zoom_extent === nothing ? 2 : 3
    # Grid position of the profile/interior/zoom axes for phase index i, and
    # of the shared colorbar/title, transposed between the two orientations.
    profile_pos(i) = orientation == :vertical ? (i, 1) : (1, i)
    quiver_pos(i)  = orientation == :vertical ? (i, 2) : (2, i)
    zoom_pos(i)    = orientation == :vertical ? (i, 3) : (3, i)
    cbar_pos       = orientation == :vertical ? (div(n, 2):(div(n, 2) + 1), n_panels + 1) : (1:n_panels, n + 1)
    title_pos      = orientation == :vertical ? (1, 1:n_panels) : (1, 1:n)

    x_prof_s = x_prof ./ length_scale
    x_grid_s = x_grid ./ length_scale
    z_grid_s = z_grid ./ length_scale
    xs_s     = xs[keep] ./ length_scale
    zs_s     = zs[keep] ./ length_scale

    dx_s      = x_grid_s[2] - x_grid_s[1]
    margin_s  = z_grid_s[1]
    max_len_s = min(2.2 * dx_s, 1.7 * margin_s)

    # Fixed-length direction-only arrows for the :contour panels: magnitude
    # is already the fill color, so these just mark where the velocity
    # points, skipping any point where it's exactly zero (direction
    # undefined -- happens right at the antisymmetric profile's zero
    # crossing, e.g. x=0 at depth). Single flat colour, magenta: :turbo runs
    # dark purple -> blue -> cyan -> green -> yellow -> orange -> dark red,
    # never through magenta/pink, so this hue reads clearly against every
    # part of the colormap rather than only the dark or only the light end.
    function draw_direction_arrows!(ax, xs_pts, zs_pts, us_pts, ws_pts, arrow_len)
        spd = sqrt.(us_pts .^ 2 .+ ws_pts .^ 2)
        valid = spd .> 0
        arrows!(ax, xs_pts[valid], zs_pts[valid],
            (us_pts[valid] ./ spd[valid]) .* arrow_len,
            (ws_pts[valid] ./ spd[valid]) .* arrow_len,
            arrowsize = 16, linewidth = 3.5, lengthscale = 1.0, color = :magenta)
    end

    # Surface velocity profile V(x) for each phase, nondimensionalized.
    all_v_prof = [make_V(tau_i).(x_prof) .* vel_nondim for tau_i in taus]

    # Surface displacement profile u(x), nondimensionalized by this block's
    # own displacement scale (b for STEP, aT for CSS).
    all_u_prof = [make_U(tau_i).(x_prof) .* disp_nondim for tau_i in taus]
    log_progress("surface profiles done")

    if interior_panel == :quiver
        # Interior velocity field (vx, vz) on the (x_grid, z_grid) mesh, one
        # (us, ws) pair per phase, via the halfspace convolution (origin point
        # excluded).
        all_us = Vector{Vector{Float64}}(undef, n)
        all_ws = Vector{Vector{Float64}}(undef, n)
        for (i, tau_i) in enumerate(taus)
            V = make_V(tau_i)
            all_us[i] = (halfspace_vx.(xs, zs, Ref(V)) .* vel_nondim)[keep]
            all_ws[i] = (halfspace_vz.(xs, zs, Ref(V)) .* vel_nondim)[keep]
        end
        speeds = [sqrt.(all_us[i] .^ 2 .+ all_ws[i] .^ 2) for i in 1:n]
        # Arrow length and color both use THIS block's own range (post-origin-
        # removal, passed in as color_min/color_max) so STEP and CSS are each
        # well populated and well contrasted visually.
        len_scale_arrow = max_len_s / maximum(maximum.(speeds))
    else  # :contour
        x_grid_dense_s = x_grid_dense ./ length_scale
        z_grid_dense_s = z_grid_dense ./ length_scale

        # |v| on the dense grid, every point included (no origin exclusion):
        # contourf's fixed level range (built from color_min/color_max below)
        # clamps the near-dike spike into the top band, the same way
        # out-of-range colors are already clamped for the quiver arrows.
        speed_grids = Vector{Matrix{Float64}}(undef, n)
        for (i, tau_i) in enumerate(taus)
            speed_grids[i] = grid_or_cache("dense_$i") do
                speed_grid(x_grid_dense, z_grid_dense, make_V(tau_i); tag = "dense grid, phase $i/$n")
            end
            log_progress("dense grid, phase $i/$n done")
        end

        # Sparse-grid (x_grid, z_grid) velocity components, for direction-
        # only arrows overlaid on the filled contour -- reuses the quiver
        # grid/mask, but only the (us, ws) *direction* is drawn (fixed arrow
        # length), since magnitude is already shown by the fill color.
        all_us_dir = Vector{Vector{Float64}}(undef, n)
        all_ws_dir = Vector{Vector{Float64}}(undef, n)
        for (i, tau_i) in enumerate(taus)
            V = make_V(tau_i)
            all_us_dir[i] = (hvx.(xs, zs, Ref(V)) .* vel_nondim)[keep]
            all_ws_dir[i] = (hvz.(xs, zs, Ref(V)) .* vel_nondim)[keep]
            log_progress("direction arrows, phase $i/$n done")
        end

        if zoom_extent !== nothing
            # Same point count as the main dense grid (nx_dense x nz_dense),
            # but resampled over the much smaller zoom_extent = (x_max_m,
            # z_max_m) box -- same resolution budget, far finer sampling per
            # unit area.
            x_max_zoom_m, z_max_zoom_m = zoom_extent
            x_grid_zoom = collect(range(-x_max_zoom_m, x_max_zoom_m, length = nx_dense))
            z_grid_zoom = collect(range(0.0, z_max_zoom_m, length = nz_dense))
            x_grid_zoom_s = x_grid_zoom ./ length_scale
            z_grid_zoom_s = z_grid_zoom ./ length_scale

            speed_grids_zoom = Vector{Matrix{Float64}}(undef, n)
            for (i, tau_i) in enumerate(taus)
                speed_grids_zoom[i] = grid_or_cache("zoom_$i") do
                    speed_grid(x_grid_zoom, z_grid_zoom, make_V(tau_i); tag = "zoom grid, phase $i/$n")
                end
                log_progress("zoom grid, phase $i/$n done")
            end

            # Sparse direction-arrow grid over the same zoom box (much
            # coarser than the dense contour grid above -- just enough
            # points to show the flow pattern). nx even keeps x=0 (the one
            # truly singular column) off the sampled grid.
            nx_zoom_arrows, nz_zoom_arrows = 10, 8
            x_grid_zoom_arrows = collect(range(-x_max_zoom_m, x_max_zoom_m, length = nx_zoom_arrows))
            z_grid_zoom_arrows = collect(range(0.0, z_max_zoom_m, length = nz_zoom_arrows))
            xs_zoom_arrows = Float64[]; zs_zoom_arrows = Float64[]
            for zi in z_grid_zoom_arrows, xi in x_grid_zoom_arrows
                push!(xs_zoom_arrows, xi); push!(zs_zoom_arrows, zi)
            end
            xs_zoom_arrows_s = xs_zoom_arrows ./ length_scale
            zs_zoom_arrows_s = zs_zoom_arrows ./ length_scale

            all_us_dir_zoom = Vector{Vector{Float64}}(undef, n)
            all_ws_dir_zoom = Vector{Vector{Float64}}(undef, n)
            for (i, tau_i) in enumerate(taus)
                V = make_V(tau_i)
                all_us_dir_zoom[i] = hvx.(xs_zoom_arrows, zs_zoom_arrows, Ref(V)) .* vel_nondim
                all_ws_dir_zoom[i] = hvz.(xs_zoom_arrows, zs_zoom_arrows, Ref(V)) .* vel_nondim
                log_progress("zoom arrows, phase $i/$n done")
            end
            dx_zoom_arrows_s = (x_grid_zoom_arrows[2] - x_grid_zoom_arrows[1]) / length_scale
        end

        # Filled-contour level edges spanning [color_min, color_max]: log-
        # spaced under a log10 colorscale (matching the Colorbar's log axis
        # below), linearly spaced otherwise. More fill bands than line levels
        # -- 30 bands give smooth color gradation within each decade (e.g.
        # 0.1-1.0), but drawing all 30 as black contour lines would be too
        # dense to read, so the lines use a coarser subset.
        n_fill_levels = 70
        n_line_levels = 10
        level_edges(n) = colorscale === log10 ?
            10 .^ range(log10(color_min), log10(color_max), length = n + 1) :
            range(color_min, color_max, length = n + 1)
        fill_levels = level_edges(n_fill_levels)
        line_levels = level_edges(n_line_levels)
    end

    for (i, (tau_i, lab)) in enumerate(zip(taus, labels))
        L = guide_length_fn(tau_i)
        v_prof_nondim = all_v_prof[i]
        u_prof_nondim = all_u_prof[i]

        # Surface velocity profile v_x(x) for this phase, with the
        # characteristic length L = a*tau marked as a vertical guide.
        # "(dimensionless)" is spelled out only on the first phase -- it
        # applies to every phase, so stating it once is enough.
        ax1 = Axis(fig[profile_pos(i)...],
            xlabel = "x / $length_symbol", ylabel = "v_x / $vel_symbol",
            title  = guide_formula === nothing ? "$time_symbol = $lab" :
                "$time_symbol = $lab   (L/$length_symbol = $guide_formula = $(round(L / length_scale, sigdigits = 3)))",
            limits = ((-x_span_m / length_scale, x_span_m / length_scale), (-vel_lim, vel_lim)))
        hlines!(ax1, [0.0], color = (:gray, 0.4), linestyle = :dot)
        hlines!(ax1, [1.0], color = (:gray, 0.4), linestyle = :dot)
        hlines!(ax1, [-1.0], color = (:gray, 0.4), linestyle = :dot)
        vlines!(ax1, [0.0], color = (:gray, 0.4), linestyle = :dot)
        vlines!(ax1, [L / length_scale], color = (:firebrick, 0.6), linestyle = :dash, linewidth = 1.5)
        l_v = lines!(ax1, x_prof_s, v_prof_nondim, color = :forestgreen, linewidth = 4.5)

        # Displacement u(x) shares the panel via a twin axis on the right,
        # linked in x but with its own (symmetric) y-scale since u and v have
        # unrelated magnitudes and units.
        ax1b = Axis(fig[profile_pos(i)...],
            ylabel = "u_x / $disp_symbol",
            yaxisposition = :right,
            limits = ((-x_span_m / length_scale, x_span_m / length_scale), (-u_lim, u_lim)))
        hidespines!(ax1b)
        hidexdecorations!(ax1b)
        linkxaxes!(ax1, ax1b)
        l_u = lines!(ax1b, x_prof_s, u_prof_nondim, color = :darkorange, linewidth = 4.5, linestyle = :dash)
        if i == 1
            axislegend(ax1, [l_v, l_u], ["v_x / $vel_symbol", "u_x / $disp_symbol"],
                position = :lt, framevisible = false, labelsize = legend_fontsize)
        end

        if interior_panel == :quiver
            us, ws = all_us[i], all_ws[i]
            speed  = speeds[i]
            println("τ = $lab:  |v| / $vel_symbol range = [$(minimum(speed)), $(maximum(speed))]")

            # (x, depth) panel for the interior velocity quiver.
            ax2 = Axis(fig[quiver_pos(i)...],
                xlabel = "x / $length_symbol", ylabel = "z / $length_symbol",
                yreversed = true, aspect = DataAspect(),
                limits = ((x_grid_s[1], x_grid_s[end]), (0.0, z_grid_s[end])))

            arrows!(ax2, xs_s, zs_s, us .* len_scale_arrow, ws .* len_scale_arrow,
                arrowsize   = 20,
                linewidth   = 5.5,
                lengthscale = 1.0,
                color       = speed, colorrange = (color_min, color_max), colormap = :turbo,
                colorscale  = colorscale)

            # Reference arrow, shown once (first phase only): calibrates arrow
            # length against a labelled |v| value, using the same len_scale_arrow
            # conversion as the data arrows above, in the empty margin between
            # the surface (z=0) and the shallowest sampled row.
            if i == 1
                # Placed in axis-relative (0-1) coordinates, not data coordinates:
                # the z<z_grid_s[1] margin's DATA size varies hugely between the
                # STEP (length_scale=H) and CSS (length_scale=aT) blocks, but a
                # fixed-pixel arrow needs a fixed on-screen fraction of room
                # regardless, so relative space is the only placement that is
                # robust to both. y=1 is the top of the panel on screen (z=0,
                # the surface) regardless of the z-axis being reversed.
                arrows!(ax2, [0.70], [0.15], [0.12], [0.0],
                    arrowsize = 20, linewidth = 5.5, lengthscale = 1.0, color = :black, space = :relative)
                text!(ax2, 0.60, 0.07, text = "|v| / $vel_symbol = $ref_speed",
                    fontsize = ref_fontsize, align = (:left, :center), space = :relative)
            end
        else  # :contour
            speed_grid = speed_grids[i]
            println("τ = $lab:  |v| / $vel_symbol range = [$(minimum(speed_grid)), $(maximum(speed_grid))]")

            # (x, depth) panel for the interior |v| filled contour map.
            ax2 = Axis(fig[quiver_pos(i)...],
                xlabel = "x / $length_symbol", ylabel = "z / $length_symbol",
                yreversed = true, aspect = DataAspect(),
                limits = ((x_grid_s[1], x_grid_s[end]), (0.0, z_grid_s[end])))

            contourf!(ax2, x_grid_dense_s, z_grid_dense_s, speed_grid,
                levels = fill_levels, colormap = :turbo, colorscale = colorscale,
                extendlow = :auto, extendhigh = :auto)
            contour!(ax2, x_grid_dense_s, z_grid_dense_s, speed_grid,
                levels = line_levels, color = (:black, 0.3), linewidth = 1.0)
            draw_direction_arrows!(ax2, xs_s, zs_s, all_us_dir[i], all_ws_dir[i], 0.6 * max_len_s)

            if zoom_extent !== nothing
                # Zoomed-in panel over the near-dike peak, same fill/line
                # levels and colormap as the main contour panel (col 2) so
                # colors read identically across both -- only the sampled
                # extent differs.
                speed_grid_zoom = speed_grids_zoom[i]
                ax3 = Axis(fig[zoom_pos(i)...],
                    xlabel = "x / $length_symbol", ylabel = "z / $length_symbol",
                    yreversed = true, aspect = DataAspect(),
                    limits = ((x_grid_zoom_s[1], x_grid_zoom_s[end]), (0.0, z_grid_zoom_s[end])))

                contourf!(ax3, x_grid_zoom_s, z_grid_zoom_s, speed_grid_zoom,
                    levels = fill_levels, colormap = :turbo, colorscale = colorscale,
                    extendlow = :auto, extendhigh = :auto)
                contour!(ax3, x_grid_zoom_s, z_grid_zoom_s, speed_grid_zoom,
                    levels = line_levels, color = (:black, 0.3), linewidth = 1.0)
                draw_direction_arrows!(ax3, xs_zoom_arrows_s, zs_zoom_arrows_s,
                    all_us_dir_zoom[i], all_ws_dir_zoom[i], 0.6 * dx_zoom_arrows_s)
            end
        end
    end

    # Log-scaled color mapping needs tick labels printed as plain decimals,
    # not powers of ten; a linear mapping needs neither scale nor tickformat.
    if colorscale === log10
        # Default log-scale ticks land only on decade powers of ten (e.g.
        # ..., 0.1, 1.0), so color_max (1.5, not a power of ten) is never
        # itself labelled even though the bar's top IS color_max -- add it
        # explicitly so the colorbar's own ceiling is always shown.
        decade_ticks = 10.0 .^ (ceil(Int, log10(color_min)):floor(Int, log10(color_max)))
        cbar_ticks = sort(unique(vcat(decade_ticks, color_max)))
        Colorbar(fig[cbar_pos...], limits = (color_min, color_max), colormap = :turbo,
            scale = log10, ticks = cbar_ticks,
            tickformat = vs -> [string(round(v, sigdigits = 2)) for v in vs],
            label = colorbar_label)
    else
        Colorbar(fig[cbar_pos...], limits = (color_min, color_max), colormap = :turbo,
            label = colorbar_label)
    end
    Label(fig[title_pos..., Top()], title, fontsize = title_fontsize, padding = (0, 0, 30, 0))

    return color_min, color_max
end
