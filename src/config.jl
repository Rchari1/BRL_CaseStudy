# Configuration: defaults, YAML loading/validation, and the typed `SimParams`
# built from it. Units are SI throughout. "per width" means per metre of
# blanket width; the model is a 1 m wide strip and totals are scaled by
# `blanket.blanket_width`.
#
# EVERY PHYSICAL DEFAULT BELOW IS A PLACEHOLDER. None of them has been measured.

"""
    default_config() -> Dict{String,Any}

Nested dictionary of every configuration key with its default (placeholder)
value. YAML files are merged on top of this; unknown keys are rejected.
"""
function default_config()
    Dict{String,Any}(
        "run_name" => "baseline",
        "blanket" => Dict{String,Any}(
            "n_panels" => 8,                     # placeholder
            "panel_length" => 0.5,               # m, placeholder
            "nodes_per_panel" => 20,             # segments per panel (numerics; 25 mm converges baseline peak to ~1%)
            "areal_density" => 1.0,              # kg/m^2, placeholder
            "EA" => 1.25e5,                      # N/m per width, placeholder
            "panel_EI" => 1.0e-2,                # N*m per width, placeholder
            "blanket_width" => 1.0,              # m, scale factor for totals
            "min_cell_bend_radius" => 0.5,       # m, placeholder allowable
            "axial_damping_ratio" => 0.02,       # placeholder
            "panel_bending_damping_ratio" => 0.02, # placeholder
        ),
        "hinge" => Dict{String,Any}(
            "law" => "linear",                   # linear | tabulated
            "k_hinge" => 2.6e-2,                 # N*m/rad per width, placeholder
            "theta_rest" => 0.3pi,               # rad, placeholder fold memory
            "damping_ratio" => 0.05,             # placeholder
            "table_file" => "",                  # CSV (fold_angle_rad, moment_Nm_per_m)
        ),
        "attachments" => Dict{String,Any}(
            "k_att_root" => 1.0e5,               # N/m per width, placeholder
            "c_att_root" => 50.0,                # N*s/m per width, placeholder
            "k_att_tip" => 1.0e5,                # N/m per width, placeholder
            "c_att_tip" => 50.0,                 # N*s/m per width, placeholder
        ),
        "deployment" => Dict{String,Any}(
            "profile" => "smooth",               # smooth | trapezoid | constant_then_stop
            "T_deploy" => 30.0,                  # s, placeholder
            "accel_fraction" => 0.2,             # trapezoid: fraction of T spent accelerating (and decelerating)
            "start_ramp_fraction" => 0.1,        # constant_then_stop: smooth start ramp fraction
            "preload_strain" => 1.0e-4,          # placeholder
            "stow_gap" => 1.0e-3,                # m, placeholder
            "T_settle" => 3.0,                   # s, settle phase with the tip held (0 disables)
            "settle_mass_damping" => 4.0,        # 1/s, extra damping applied only while settling
            "T_hold" => 5.0,                     # s, simulated after the tip stops
        ),
        "environment" => Dict{String,Any}(
            "gravity" => false,
            "gravity_direction" => "-z",         # -z | +z | -x | +x
            "g" => 9.80665,                      # m/s^2
            "drag" => false,
            "air_density" => 1.2,                # kg/m^3
            "drag_coefficient" => 1.5,           # flat-plate estimate, placeholder
        ),
        "numerics" => Dict{String,Any}(
            "dt" => "auto",                      # "auto" or a fixed step in s
            "dt_safety" => 0.5,                  # fraction of the estimated critical step
            "axial_stiffness_scale" => 1.0,      # <1 softens EA for exploration; LOWERS snap peaks
        ),
        "output" => Dict{String,Any}(
            "dir" => "results",
            "field_dt" => 0.01,                  # s, shape/tension/angle snapshots
            "channel_dt" => 5.0e-4,              # s, forces/energies/peaks time series
            "snap_window_pre" => 0.03,           # s of high-rate tension field kept before the peak
            "snap_window_post" => 0.09,          # s kept after the peak
            "animation" => true,
            "animation_fps" => 30,
            "animation_seconds" => 16.0,         # video length (deployment is time-compressed)
            "gif" => true,
        ),
    )
end

"""Deep-merge `src` into `dst`, erroring on keys that `dst` does not define."""
function merge_config!(dst::Dict{String,Any}, src::AbstractDict; path = "")
    for (k, v) in src
        key = string(k)
        full = isempty(path) ? key : "$path.$key"
        haskey(dst, key) || error("Unknown configuration key `$full`")
        if dst[key] isa Dict
            v isa AbstractDict || error("Configuration key `$full` must be a mapping")
            merge_config!(dst[key], v; path = full)
        else
            dst[key] = v
        end
    end
    dst
end

"""
    load_config(path) -> Dict{String,Any}

Read a YAML file and merge it over [`default_config`](@ref). A key
`base: other.yaml` (relative to the file) is merged first, so configs can
inherit from the baseline.
"""
function load_config(path::AbstractString)
    raw = YAML.load_file(path; dicttype = Dict{String,Any})
    raw === nothing && (raw = Dict{String,Any}())
    cfg = if haskey(raw, "base")
        load_config(joinpath(dirname(path), pop!(raw, "base")))
    else
        default_config()
    end
    haskey(raw, "run_name") || (cfg["run_name"] = splitext(basename(path))[1])
    merge_config!(cfg, raw)
end

"""Read a dotted key such as `"hinge.k_hinge"`."""
function get_param(cfg::Dict, dotted::AbstractString)
    node = cfg
    for part in split(dotted, '.')
        haskey(node, part) || error("Unknown configuration key `$dotted`")
        node = node[part]
    end
    node
end

"""Set a dotted key such as `"hinge.k_hinge"` (the key must already exist)."""
function set_param!(cfg::Dict, dotted::AbstractString, value)
    parts = split(dotted, '.')
    node = cfg
    for part in parts[1:end-1]
        haskey(node, part) || error("Unknown configuration key `$dotted`")
        node = node[part]
    end
    haskey(node, parts[end]) || error("Unknown configuration key `$dotted`")
    node[parts[end]] = value
    cfg
end

"""Parse numbers or strings like `"0.3pi"`, `"0.3*pi"`, `"pi"`, `"1e-3"`."""
function parse_number(v, key)
    v isa Real && return Float64(v)
    if v isa AbstractString
        s = replace(lowercase(strip(v)), " " => "", "π" => "pi")
        m = match(r"^([-+0-9.e]*)\*?pi$", s)
        if m !== nothing
            c = m.captures[1]
            coeff = c in ("", "+") ? 1.0 : c == "-" ? -1.0 : parse(Float64, c)
            return coeff * pi
        end
        x = tryparse(Float64, s)
        x === nothing || return x
    end
    error("Configuration key `$key` must be a number (got $(repr(v)))")
end

num(cfg, key) = parse_number(get_param(cfg, key), key)

function positive(cfg, key; allow_zero = false)
    x = num(cfg, key)
    ok = allow_zero ? x >= 0 : x > 0
    ok || error("Configuration key `$key` must be $(allow_zero ? "non-negative" : "positive") (got $x)")
    x
end

const GRAVITY_DIRECTIONS = Dict(
    "-z" => SVector(0.0, -1.0), "+z" => SVector(0.0, 1.0),
    "-x" => SVector(-1.0, 0.0), "+x" => SVector(1.0, 0.0),
)

"""
    SimParams

Typed, validated, derived parameters for one run. Built from a config
dictionary by `SimParams(cfg)`. Stiffness and damping values are per metre
of blanket width; masses are kg per metre of width.
"""
struct SimParams{H<:HingeLaw,P<:DeployProfile}
    cfg::Dict{String,Any}
    name::String
    layout::Layout
    # material (per width)
    EA::Float64                  # effective EA after axial_stiffness_scale
    axial_stiffness_scale::Float64
    k_axial::Vector{Float64}     # per segment spring, N/m
    c_axial::Vector{Float64}     # per segment dashpot, N*s/m
    k_bend::Vector{Float64}      # per node torsional stiffness (0 at ends and hinges)
    c_bend::Vector{Float64}      # per node torsional damping (panel interior)
    hinge_law::H
    theta_rest::Float64          # fold-angle rest value used by the linear law
    fold_sign::Vector{Float64}   # per node: +-1 at hinges (direction of the stowed fold), 0 elsewhere
    min_cell_bend_radius::Float64
    width::Float64
    # attachments
    k_root::Float64
    c_root::Float64
    k_tip::Float64
    c_tip::Float64
    # deployment
    profile::P
    preload_strain::Float64
    stow_gap::Float64
    T_settle::Float64
    T_deploy::Float64
    T_hold::Float64
    settle_alpha::Float64
    # environment
    gravity::SVector{2,Float64}  # zero vector when disabled
    drag::Bool
    air_density::Float64
    drag_coefficient::Float64
    # numerics
    dt::Float64
    dt_crit::Float64
    # output
    field_dt::Float64
    channel_dt::Float64
    snap_pre::Float64
    snap_post::Float64
end

function SimParams(cfg::Dict{String,Any})
    n_panels = Int(num(cfg, "blanket.n_panels"))
    n_panels >= 1 || error("blanket.n_panels must be >= 1")
    nps = Int(num(cfg, "blanket.nodes_per_panel"))
    nps >= 2 || error("blanket.nodes_per_panel must be >= 2")
    Lp = positive(cfg, "blanket.panel_length")
    rho = positive(cfg, "blanket.areal_density")
    scale = positive(cfg, "numerics.axial_stiffness_scale")
    EA = positive(cfg, "blanket.EA") * scale
    EI = positive(cfg, "blanket.panel_EI")
    lay = build_layout(n_panels, Lp, nps, rho)
    L_total = n_panels * Lp

    # Axial: stiffness-proportional (material) damping, damping ratio defined
    # at the fundamental axial mode of the fully deployed blanket (fixed-fixed).
    zeta_ax = positive(cfg, "blanket.axial_damping_ratio"; allow_zero = true)
    omega_ax = pi * sqrt(EA / rho) / L_total
    eta_ax = 2zeta_ax * EA / omega_ax                  # N*s per width
    k_axial = EA ./ lay.rest_length
    c_axial = eta_ax ./ lay.rest_length

    # Panel bending: discrete EI/l_bar at interior nodes; stiffness-proportional
    # damping with the ratio defined at the first free-free panel bending mode.
    zeta_b = positive(cfg, "blanket.panel_bending_damping_ratio"; allow_zero = true)
    omega_b = 22.373 * sqrt(EI / (rho * Lp^4))
    k_bend = zeros(lay.N)
    c_bend = zeros(lay.N)
    for i in 2:lay.N-1
        lay.kind[i] == PANEL_NODE || continue
        lbar = 0.5 * (lay.rest_length[i-1] + lay.rest_length[i])
        k_bend[i] = EI / lbar
        c_bend[i] = 2zeta_b * k_bend[i] / omega_b
    end

    # Hinges: damping ratio defined against one panel rotating about the hinge.
    k_h = positive(cfg, "hinge.k_hinge"; allow_zero = true)
    theta_rest = num(cfg, "hinge.theta_rest")
    zeta_h = positive(cfg, "hinge.damping_ratio"; allow_zero = true)
    I_panel = rho * Lp^3 / 3
    c_h = 2zeta_h * sqrt(k_h * I_panel)
    law_name = lowercase(string(get_param(cfg, "hinge.law")))
    hinge_law = if law_name == "linear"
        LinearHinge(k_h, theta_rest, c_h)
    elseif law_name == "tabulated"
        file = string(get_param(cfg, "hinge.table_file"))
        isempty(file) && error("hinge.law = tabulated requires hinge.table_file")
        TabulatedHinge(file, c_h)
    else
        error("hinge.law must be `linear` or `tabulated` (got `$law_name`)")
    end

    fold_sign = zeros(lay.N)
    for (k, i) in enumerate(lay.hinge_nodes)
        fold_sign[i] = isodd(k) ? 1.0 : -1.0   # matches stowed_positions()
    end

    profile = make_profile(string(get_param(cfg, "deployment.profile")),
        positive(cfg, "deployment.T_deploy"),
        num(cfg, "deployment.accel_fraction"),
        num(cfg, "deployment.start_ramp_fraction"))

    gravity = if get_param(cfg, "environment.gravity") == true
        dir = string(get_param(cfg, "environment.gravity_direction"))
        haskey(GRAVITY_DIRECTIONS, dir) || error("environment.gravity_direction must be one of $(sort(collect(keys(GRAVITY_DIRECTIONS))))")
        positive(cfg, "environment.g") * GRAVITY_DIRECTIONS[dir]
    else
        SVector(0.0, 0.0)
    end

    k_root = positive(cfg, "attachments.k_att_root"; allow_zero = true)
    c_root = positive(cfg, "attachments.c_att_root"; allow_zero = true)
    k_tip = positive(cfg, "attachments.k_att_tip"; allow_zero = true)
    c_tip = positive(cfg, "attachments.c_att_tip"; allow_zero = true)
    settle_alpha = positive(cfg, "deployment.settle_mass_damping"; allow_zero = true)

    dt_crit = critical_dt(lay, k_axial, c_axial, k_bend, c_bend, hinge_law,
        k_root, c_root, k_tip, c_tip, settle_alpha)
    dt_cfg = get_param(cfg, "numerics.dt")
    dt = if dt_cfg isa AbstractString && lowercase(dt_cfg) == "auto"
        positive(cfg, "numerics.dt_safety") * dt_crit
    else
        positive(cfg, "numerics.dt")
    end

    SimParams(cfg, string(get_param(cfg, "run_name")), lay, EA, scale, k_axial, c_axial,
        k_bend, c_bend, hinge_law, theta_rest, fold_sign,
        positive(cfg, "blanket.min_cell_bend_radius"), positive(cfg, "blanket.blanket_width"),
        k_root, c_root, k_tip, c_tip, profile,
        positive(cfg, "deployment.preload_strain"; allow_zero = true),
        positive(cfg, "deployment.stow_gap"),
        positive(cfg, "deployment.T_settle"; allow_zero = true),
        positive(cfg, "deployment.T_deploy"),
        positive(cfg, "deployment.T_hold"; allow_zero = true),
        settle_alpha, gravity, get_param(cfg, "environment.drag") == true,
        positive(cfg, "environment.air_density"; allow_zero = true),
        positive(cfg, "environment.drag_coefficient"; allow_zero = true),
        dt, dt_crit,
        positive(cfg, "output.field_dt"), positive(cfg, "output.channel_dt"),
        positive(cfg, "output.snap_window_pre"; allow_zero = true),
        positive(cfg, "output.snap_window_post"; allow_zero = true))
end

SimParams(path::AbstractString) = SimParams(load_config(path))

"""Total flat (unstretched) blanket length, m."""
flat_length(p::SimParams) = p.layout.n_panels * p.layout.panel_length

"""Blanket mass per metre of width, kg/m."""
blanket_mass(p::SimParams) = sum(p.layout.mass)

"""
    print_startup(io, p)

Print the stability estimate, wave speed and any warnings about settings that
bias the results.
"""
function print_startup(io::IO, p::SimParams)
    lay = p.layout
    c = sqrt(p.EA / lay.areal_density)
    @printf(io, "[%s] %d panels x %.3g m, %d nodes, %s profile, T_deploy = %.3g s\n",
        p.name, lay.n_panels, lay.panel_length, lay.N, profile_name(p.profile), p.T_deploy)
    @printf(io, "  axial wave speed sqrt(EA/rho) = %.1f m/s; critical dt ~ %.3e s; using dt = %.3e s (%.2f x critical)\n",
        c, p.dt_crit, p.dt, p.dt / p.dt_crit)
    if p.axial_stiffness_scale != 1.0
        @printf(io, "  WARNING: axial_stiffness_scale = %.3g (EA softened). Snap-taut peak loads are UNDER-predicted.\n",
            p.axial_stiffness_scale)
    end
    if p.dt > p.dt_crit
        println(io, "  WARNING: dt exceeds the stability estimate; the run may diverge.")
    end
end
