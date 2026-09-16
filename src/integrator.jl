# Time integration: symplectic (semi-implicit) Euler, i.e. leapfrog with
# velocities at half steps, fixed step from a stability estimate.
#
#   v_{n+1/2} = v_{n-1/2} + dt * F(x_n, v_{n-1/2}, t_n) / m
#   x_{n+1}   = x_n + dt * v_{n+1/2}
#
# Energy bookkeeping uses the exact leapfrog identity
#   1/2 m|v_{n+1/2}|^2 - 1/2 m|v_{n-1/2}|^2 = dt F_n . (v_{n-1/2} + v_{n+1/2})/2,
# so dissipation and attachment work are integrated with the same mid-step
# velocity the scheme itself uses. The reported balance residual then measures
# only the O(dt^2) error in the conservative potentials.

"""
    critical_dt(...)

Gershgorin-type estimate of the explicit stability limit. For every node it
bounds the highest local frequency from the row sums of the stiffness
(axial, bending incl. a geometric term, attachment) and damping matrices, and
applies the damped central-difference limit dt < 2/w (sqrt(1+z^2) - z).
"""
function critical_dt(lay::Layout, k_axial, c_axial, k_bend, c_bend, law::HingeLaw,
    k_root, c_root, k_tip, c_tip, settle_alpha)
    N = lay.N
    K = zeros(N)
    C = zeros(N)
    for j in 1:N-1
        K[j] += 2k_axial[j]; K[j+1] += 2k_axial[j]
        C[j] += 2c_axial[j]; C[j+1] += 2c_axial[j]
    end
    kh = hinge_stiffness_bound(law)
    ch = hinge_damping(law)
    for i in 2:N-1
        la, lb = lay.rest_length[i-1], lay.rest_length[i]
        g = (1 / la, 1 / la + 1 / lb, 1 / lb)       # |dtheta/dx| at i-1, i, i+1
        G = sum(g)
        if lay.kind[i] == PANEL_NODE
            kk, cc, geo = k_bend[i], c_bend[i], 0.0
        else
            kk, cc, geo = kh, ch, kh * pi            # geometric term from a moment up to k*pi
        end
        for (off, gi) in zip((-1, 0, 1), g)
            K[i+off] += kk * gi * G + geo * 2gi * G / 3
            C[i+off] += cc * gi * G
        end
    end
    K[1] += k_root; C[1] += c_root
    K[N] += k_tip; C[N] += c_tip
    dt = Inf
    for i in 1:N
        m = lay.mass[i]
        w = sqrt(K[i] / m)
        gamma = C[i] / m + settle_alpha
        z = gamma / (2w)
        dt = min(dt, 2 / w * (sqrt(1 + z^2) - z), 2 / gamma)
    end
    dt
end

"""Scratch buffers for one run (not shared between threads)."""
struct Workspace
    F::Vector{SVector{2,Float64}}
    Fd_axial::Vector{SVector{2,Float64}}
    Fd_panel::Vector{SVector{2,Float64}}
    Fd_hinge::Vector{SVector{2,Float64}}
    Fd_drag::Vector{SVector{2,Float64}}
    Fd_settle::Vector{SVector{2,Float64}}
    vnew::Vector{SVector{2,Float64}}
    vbar::Vector{SVector{2,Float64}}
    theta::Vector{Float64}
    curvature::Vector{Float64}
    tension::Vector{Float64}
end

function Workspace(N::Int)
    z() = zeros(SVector{2,Float64}, N)
    Workspace(z(), z(), z(), z(), z(), z(), z(), z(), zeros(N), zeros(N), zeros(N - 1))
end

"""Sum of F[i] . v[i] (power of a nodal force set)."""
function power(F, v)
    s = 0.0
    @inbounds @simd for i in eachindex(F)
        s += dot(F[i], v[i])
    end
    s
end

"""Index of the maximum, the maximum and the minimum of a vector in one pass."""
function extremes(x)
    jmax, hi, lo = 1, x[1], x[1]
    @inbounds for j in 2:length(x)
        xj = x[j]
        if xj > hi
            hi, jmax = xj, j
        end
        lo = min(lo, xj)
    end
    (jmax, hi, lo)
end

"""Energies and attachment forces produced by one force evaluation."""
struct ForceEval
    E_axial::Float64
    E_panel::Float64
    E_hinge::Float64
    E_att::Float64
    PE_grav::Float64
    Fs_root::SVector{2,Float64}   # forces ON the root node (N/m)
    Fd_root::SVector{2,Float64}
    Fs_tip::SVector{2,Float64}
    Fd_tip::SVector{2,Float64}
end

"""
    evaluate_forces!(ws, p, x, v, t, root_anchor, tip_anchor, tip_velocity) -> ForceEval

Assemble all nodal forces at state (x, v) and time t (s; negative while settling).
"""
function evaluate_forces!(ws::Workspace, p::SimParams, x, v, t::Float64,
    root_anchor::SVector{2,Float64}, a_tip::SVector{2,Float64}, adot_tip::SVector{2,Float64})
    lay = p.layout
    N = lay.N
    zero2 = SVector(0.0, 0.0)
    fill!(ws.F, zero2)
    fill!(ws.Fd_axial, zero2); fill!(ws.Fd_panel, zero2); fill!(ws.Fd_hinge, zero2)
    p.drag && fill!(ws.Fd_drag, zero2)
    t < 0 && fill!(ws.Fd_settle, zero2)

    E_axial = axial_forces!(ws.F, ws.Fd_axial, ws.tension, x, v, lay.rest_length, p.k_axial, p.c_axial)
    E_panel, E_hinge = bending_forces!(ws.F, ws.Fd_panel, ws.Fd_hinge, ws.theta, ws.curvature,
        x, v, lay, p.k_bend, p.c_bend, p.hinge_law, p.fold_sign)

    Fs_root, Fd_root = attachment_force(x[1], v[1], root_anchor, zero2, p.k_root, p.c_root)
    Fs_tip, Fd_tip = attachment_force(x[N], v[N], a_tip, adot_tip, p.k_tip, p.c_tip)
    ws.F[1] += Fs_root + Fd_root
    ws.F[N] += Fs_tip + Fd_tip
    E_att = attachment_energy(x[1], root_anchor, p.k_root) + attachment_energy(x[N], a_tip, p.k_tip)

    PE = 0.0
    if p.gravity != zero2
        PE = gravity_forces!(ws.F, lay.mass, p.gravity, x)
    end
    if p.drag
        drag_forces!(ws.F, ws.Fd_drag, x, v, lay.tributary, p.air_density, p.drag_coefficient)
    end
    if t < 0 && p.settle_alpha > 0
        mass_damping_forces!(ws.F, ws.Fd_settle, lay.mass, v, p.settle_alpha)
    end
    ForceEval(E_axial, E_panel, E_hinge, E_att, PE, Fs_root, Fd_root, Fs_tip, Fd_tip)
end

# Names of the recorded time-series channels (sampled every `channel_dt`).
const CHANNELS = [
    "F_root_x", "F_root_z", "F_tip_x", "F_tip_z",   # load ON the structure, N/m width
    "F_root", "F_tip",                              # magnitudes, N/m width
    "tension_max", "tension_max_arc",               # N/m width, m from root
    "curvature_ratio_max",                          # max |kappa| * R_min over panel interiors
    "KE", "E_axial", "E_panel", "E_hinge", "E_att", "PE_grav",
    "W_tip", "D_axial", "D_panel", "D_hinge", "D_att", "D_drag", "D_settle",
    "balance_residual",                             # (E - E0) - W + D, J/m width
    "anchor_x", "anchor_z", "tip_x", "tip_z", "progress",
]

"""
    Recorder

Decimated outputs: channels every `channel_dt`; shape, tension, fold-angle and
curvature fields every `field_dt`; and a high-rate tension window around the
peak attachment force (a ring buffer is copied when the window completes).
"""
mutable struct Recorder
    t_ch::Vector{Float64}
    ch::Matrix{Float64}            # (n_samples, n_channels)
    n_ch::Int
    next_ch::Float64
    t_f::Vector{Float64}
    pos::Array{Float64,3}          # (n_frames, N, 2)
    tension::Matrix{Float64}       # (n_frames, N-1)
    fold::Matrix{Float64}          # (n_frames, n_hinges)
    curv::Matrix{Float64}          # (n_frames, N) curvature ratio, signed
    n_f::Int
    next_f::Float64
    # snap window
    ring_t::Vector{Float64}
    ring_T::Matrix{Float64}
    ring_Froot::Vector{Float64}
    ring_Ftip::Vector{Float64}
    ring_head::Int                 # total frames written
    capture_at::Int                # frame count at which to copy the ring (0 = none)
    snap_t::Vector{Float64}
    snap_T::Matrix{Float64}
    snap_Froot::Vector{Float64}
    snap_Ftip::Vector{Float64}
    n_post::Int
end

function Recorder(p::SimParams, t0, t1; record_fields = true)
    lay = p.layout
    N = lay.N
    nch = floor(Int, (t1 - t0) / p.channel_dt) + 2
    nf = record_fields ? floor(Int, (t1 - t0) / p.field_dt) + 2 : 0
    npre = round(Int, p.snap_pre / p.channel_dt)
    npost = round(Int, p.snap_post / p.channel_dt)
    nring = npre + npost + 1
    Recorder(zeros(nch), zeros(nch, length(CHANNELS)), 0, t0,
        zeros(nf), zeros(nf, N, 2), zeros(nf, N - 1), zeros(nf, length(lay.hinge_nodes)), zeros(nf, N), 0, t0,
        zeros(nring), zeros(nring, N - 1), zeros(nring), zeros(nring), 0, 0,
        Float64[], zeros(0, N - 1), Float64[], Float64[], npost)
end

function copy_ring!(r::Recorder)
    n = length(r.ring_t)
    count = min(r.ring_head, n)
    idx = [mod1(r.ring_head - count + k, n) for k in 1:count]
    r.snap_t = r.ring_t[idx]
    r.snap_T = r.ring_T[idx, :]
    r.snap_Froot = r.ring_Froot[idx]
    r.snap_Ftip = r.ring_Ftip[idx]
end

"""Peak values tracked at every time step during deployment (t >= 0)."""
mutable struct Peaks
    F_root::Float64
    t_root::Float64
    F_tip::Float64
    t_tip::Float64
    tension::Float64
    t_tension::Float64
    arc_tension::Float64
    curvature_ratio::Float64
    t_curvature::Float64
    arc_curvature::Float64
    min_tension::Float64
end
Peaks() = Peaks(0.0, NaN, 0.0, NaN, -Inf, NaN, NaN, 0.0, NaN, NaN, Inf)

"""
    SimResult

Everything one run produces. Channels are loads/energies vs time; fields are
decimated snapshots; `snap_*` is the high-rate tension window around the
attachment-force peak.
"""
struct SimResult
    params::SimParams
    t_ch::Vector{Float64}
    ch::Dict{String,Vector{Float64}}
    t_f::Vector{Float64}
    pos::Array{Float64,3}
    tension::Matrix{Float64}
    fold::Matrix{Float64}
    curv::Matrix{Float64}
    snap_t::Vector{Float64}
    snap_T::Matrix{Float64}
    snap_Froot::Vector{Float64}
    snap_Ftip::Vector{Float64}
    peaks::Peaks
    E0::Float64
    balance_error::Float64         # max |residual| / energy scale
    steps::Int
    wall_time::Float64
    diverged::Bool
end

"""
    simulate(p; x0, tip_motion, t_start, t_end, record_fields, verbose) -> SimResult

Run one deployment. By default: stowed initial condition, a settle phase of
`T_settle` s (t < 0) with the tip held, then the prescribed deployment from
t = 0 to `T_deploy + T_hold`. Tests pass their own `x0` / `tip_motion`.
"""
function simulate(p::SimParams;
    x0 = stowed_positions(p.layout, p.stow_gap),
    v0 = nothing,
    root_anchor = x0[1],
    tip_motion = AnchorMotion(x0[end], root_anchor + deployed_tip_position(p), p.profile),
    t_start = -p.T_settle,
    t_end = p.T_deploy + p.T_hold,
    record_fields = true,
    verbose = true,
    io = stdout)

    verbose && print_startup(io, p)
    lay = p.layout
    N = lay.N
    m = lay.mass
    dt = p.dt
    x = collect(SVector{2,Float64}, x0)
    v = v0 === nothing ? zeros(SVector{2,Float64}, N) : collect(SVector{2,Float64}, v0)
    ws = Workspace(N)
    ws.theta .= turning_angles(x)
    rec = Recorder(p, t_start, t_end; record_fields)
    pk = Peaks()
    interior = panel_interior_nodes(lay)
    Rmin = p.min_cell_bend_radius

    nsteps = ceil(Int, (t_end - t_start) / dt - 1e-9)
    W = 0.0
    D = zero(SVector{6,Float64})   # axial, panel, hinge, att, drag, settle
    E0 = NaN
    max_resid = 0.0
    scale = 0.0
    diverged = false
    snap_pending = false
    peak_any = 0.0
    wall0 = time()
    t = t_start
    report_every = max(1, nsteps ÷ 10)

    for n in 0:nsteps
        t = t_start + n * dt
        a_tip, adot_tip = anchor_state(tip_motion, t)
        fe = evaluate_forces!(ws, p, x, v, t, root_anchor, a_tip, adot_tip)

        # Velocity update and energy accounting with the mid-step velocity.
        KE = 0.0
        @inbounds @simd for i in 1:N
            vn = v[i] + (dt / m[i]) * ws.F[i]
            ws.vnew[i] = vn
            ws.vbar[i] = 0.5 * (v[i] + vn)
            KE += 0.5 * m[i] * dot(v[i], vn)
        end
        vbar_root = ws.vbar[1]
        vbar_tip = ws.vbar[N]
        P = SVector(-power(ws.Fd_axial, ws.vbar), -power(ws.Fd_panel, ws.vbar), -power(ws.Fd_hinge, ws.vbar),
            -dot(fe.Fd_root, vbar_root) - dot(fe.Fd_tip, vbar_tip - adot_tip),
            p.drag ? -power(ws.Fd_drag, ws.vbar) : 0.0,
            t < 0 ? -power(ws.Fd_settle, ws.vbar) : 0.0)
        P_in = dot(fe.Fs_tip + fe.Fd_tip, adot_tip)

        # Values centred on t_n (half of this step's power is included).
        Wc = W + 0.5dt * P_in
        Dc = D + 0.5dt * P
        E = KE + fe.E_axial + fe.E_panel + fe.E_hinge + fe.E_att + fe.PE_grav
        isnan(E0) && (E0 = E)
        resid = (E - E0) - Wc + sum(Dc)
        scale = max(scale, abs(Wc), sum(Dc), abs(E - E0))
        max_resid = max(max_resid, abs(resid))
        W += dt * P_in
        D += dt * P

        # Loads on the structure (negative of force on the node).
        L_root = -(fe.Fs_root + fe.Fd_root)
        L_tip = -(fe.Fs_tip + fe.Fd_tip)
        Froot = norm(L_root)
        Ftip = norm(L_tip)
        jT, Tmax, Tmin = extremes(ws.tension)
        kmax = 0.0
        ikmax = 1
        @inbounds for i in interior
            c = abs(ws.curvature[i])
            if c > kmax
                kmax = c
                ikmax = i
            end
        end
        kmax *= Rmin

        if t >= 0
            if Froot > pk.F_root
                pk.F_root, pk.t_root = Froot, t
            end
            if Ftip > pk.F_tip
                pk.F_tip, pk.t_tip = Ftip, t
            end
            if max(Froot, Ftip) > peak_any
                peak_any = max(Froot, Ftip)
                snap_pending = true
            end
            if Tmax > pk.tension
                pk.tension, pk.t_tension = Tmax, t
                pk.arc_tension = 0.5 * (lay.arc[jT] + lay.arc[jT+1])
            end
            pk.min_tension = min(pk.min_tension, Tmin)
            if kmax > pk.curvature_ratio
                pk.curvature_ratio, pk.t_curvature, pk.arc_curvature = kmax, t, lay.arc[ikmax]
            end
        end

        # Channel sample (and snap-window ring buffer at the same rate).
        if t >= rec.next_ch - 0.5dt && rec.n_ch < length(rec.t_ch)
            rec.n_ch += 1
            k = rec.n_ch
            rec.t_ch[k] = t
            s, _ = progress(p.profile, t)
            vals = (L_root[1], L_root[2], L_tip[1], L_tip[2], Froot, Ftip,
                Tmax, 0.5 * (lay.arc[jT] + lay.arc[jT+1]), kmax,
                KE, fe.E_axial, fe.E_panel, fe.E_hinge, fe.E_att, fe.PE_grav,
                Wc, Dc[1], Dc[2], Dc[3], Dc[4], Dc[5], Dc[6], resid,
                a_tip[1], a_tip[2], x[N][1], x[N][2], s)
            for (c, val) in enumerate(vals)
                rec.ch[k, c] = val
            end
            rec.next_ch = t_start + k * p.channel_dt

            nring = length(rec.ring_t)
            rec.ring_head += 1
            slot = mod1(rec.ring_head, nring)
            rec.ring_t[slot] = t
            rec.ring_T[slot, :] .= ws.tension
            rec.ring_Froot[slot] = Froot
            rec.ring_Ftip[slot] = Ftip
            if snap_pending
                rec.capture_at = rec.ring_head + rec.n_post
                snap_pending = false
            end
            if rec.capture_at > 0 && rec.ring_head == rec.capture_at
                copy_ring!(rec)
                rec.capture_at = 0
            end
        end

        # Field snapshot.
        if record_fields && t >= rec.next_f - 0.5dt && rec.n_f < length(rec.t_f)
            rec.n_f += 1
            k = rec.n_f
            rec.t_f[k] = t
            @inbounds for i in 1:N
                rec.pos[k, i, 1] = x[i][1]
                rec.pos[k, i, 2] = x[i][2]
                rec.curv[k, i] = ws.curvature[i] * Rmin
            end
            rec.tension[k, :] .= ws.tension
            for (h, i) in enumerate(lay.hinge_nodes)
                rec.fold[k, h] = p.fold_sign[i] * ws.theta[i]
            end
            rec.next_f = t_start + k * p.field_dt
        end

        if !isfinite(KE) || KE > 1e12
            diverged = true
            verbose && println(io, "  DIVERGED at t = $(t) s")
            break
        end

        # Position update.
        @inbounds for i in 1:N
            v[i] = ws.vnew[i]
            x[i] = x[i] + dt * v[i]
        end
        if verbose && n > 0 && n % report_every == 0
            @printf(io, "  t = %7.3f s  (%3.0f%%)  tip load %.3g N/m\r", t, 100n / nsteps, Ftip)
        end
    end
    rec.capture_at > 0 && copy_ring!(rec)
    wall = time() - wall0
    verbose && @printf(io, "  %d steps in %.1f s wall (%.2f us/step)%s\n", nsteps, wall, 1e6wall / max(nsteps, 1), " "^20)

    nch = rec.n_ch
    ch = Dict(name => rec.ch[1:nch, c] for (c, name) in enumerate(CHANNELS))
    nf = rec.n_f
    SimResult(p, rec.t_ch[1:nch], ch, rec.t_f[1:nf], rec.pos[1:nf, :, :], rec.tension[1:nf, :],
        rec.fold[1:nf, :], rec.curv[1:nf, :], rec.snap_t, rec.snap_T, rec.snap_Froot, rec.snap_Ftip,
        pk, E0, max_resid / max(scale, eps()), nsteps, wall, diverged)
end
