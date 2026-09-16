# Figures and animation for one run. Every figure carries the run context and
# the placeholder-parameter disclaimer; forces are per metre of width with the
# scaled total stated alongside.

using CairoMakie
using CairoMakie: Makie
using .Style

const S = Style

# ----------------------------------------------------------------- helpers --

"""Format a number with ~3 significant figures and thousands separators."""
function fmt(x::Real; sig = 3)
    isfinite(x) || return string(x)
    ax = abs(x)
    ax < 1e-4 && return "0"
    if ax >= 1e5 || ax < 1e-3
        return @sprintf("%.2e", x)
    end
    digits = max(0, sig - 1 - floor(Int, log10(ax)))
    s = string(round(x; digits))
    if digits == 0
        s = string(round(Int, x))
    end
    if '.' in s
        s = rstrip(rstrip(s, '0'), '.')
    end
    intpart, frac = occursin('.', s) ? split(s, '.'; limit = 2) : (s, "")
    neg = startswith(intpart, "-")
    body = neg ? intpart[2:end] : intpart
    grouped = reverse(join(Iterators.partition(reverse(body), 3) .|> join, ","))
    (neg ? "-" : "") * grouped * (isempty(frac) ? "" : "." * frac)
end

"""Scientific notation with a real multiplication sign and superscripts: 1.2 × 10⁻⁶."""
function sci(x::Real; digits = 1)
    x == 0 && return "0"
    e = floor(Int, log10(abs(x)))
    m = x / 10.0^e
    sup = Dict('-' => '⁻', '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
        '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹')
    string(round(m; digits), " × 10", map(c -> sup[c], string(e)))
end

function env_label(p::SimParams)
    grav = p.gravity != SVector(0.0, 0.0)
    if !grav && !p.drag
        "in orbit (no gravity, no air)"
    else
        parts = String[]
        grav && push!(parts, "gravity " * string(get_param(p.cfg, "environment.gravity_direction")))
        p.drag && push!(parts, "air drag")
        "ground test (" * join(parts, ", ") * ")"
    end
end

run_context(p::SimParams) = @sprintf("%s  ·  %s profile, T_deploy = %s s  ·  %d × %s m panels  ·  %s",
    p.name, profile_name(p.profile), fmt(p.T_deploy), p.layout.n_panels, fmt(p.layout.panel_length), env_label(p))

const DISCLAIMER = "Placeholder-parameter estimate. No input has been measured; not validated against test data. Use for trends and scaling only."

function header!(pos, title, subtitle, context = "")
    g = GridLayout(pos; tellwidth = false)
    Label(g[1, 1], title; fontsize = 20, font = :medium, color = S.INK, halign = :left, tellwidth = false)
    Label(g[2, 1], subtitle; fontsize = 12.5, color = S.INK2, halign = :left, tellwidth = false)
    isempty(context) || Label(g[3, 1], context; fontsize = 11.5, color = S.MUTED, halign = :left, tellwidth = false)
    rowgap!(g, 3)
    g
end

function footer!(pos, p::SimParams; extra = "")
    g = GridLayout(pos; tellwidth = false)
    ax = Axis(g[1, 1]; width = 12, height = 12, backgroundcolor = :transparent)
    hidedecorations!(ax); hidespines!(ax)
    scatter!(ax, [0.0], [0.0]; marker = :utriangle, markersize = 11, color = S.WARNING)
    txt = DISCLAIMER
    p.axial_stiffness_scale != 1 && (txt *= @sprintf("  EA softened ×%s: snap peaks are under-predicted.", fmt(p.axial_stiffness_scale)))
    isempty(extra) || (txt *= "  " * extra)
    Label(g[1, 2], txt; fontsize = 11, color = S.INK2, halign = :left, tellwidth = false)
    colgap!(g, 6)
    g
end

"""Load ticks for a pseudo-log axis spanning 0..ymax."""
function log_ticks(ymax)
    vals = [0.0]
    labels = ["0"]
    first = ymax >= 10 ? 0 : ymax >= 1 ? -1 : -2
    for e in first:5
        v = 10.0^e
        v > ymax * 1.05 && break
        push!(vals, v)
        push!(labels, fmt(v))
    end
    (vals, labels)
end

"""Direct label: short colored key + ink text (text never wears the series color)."""
function keyed_label!(ax, x, y, text, color; align = (:left, :center), offset = (8, 0), fontsize = 11.5)
    scatter!(ax, [x], [y]; marker = :hline, markersize = 14, color, strokewidth = 0)
    text!(ax, x, y; text, align, offset, fontsize, color = S.INK2)
end

function peak_marker!(ax, x, y, color)
    scatter!(ax, [x], [y]; color, markersize = 10, strokecolor = S.SURFACE, strokewidth = 2)
end

"""Min/max envelope decimation that preserves peaks (for dense time series)."""
function envelope(t, y, nbins)
    n = length(t)
    n <= 2nbins && return (t, y)
    edges = round.(Int, range(1, n + 1; length = nbins + 1))
    tt = Float64[]; yy = Float64[]
    for b in 1:nbins
        i0, i1 = edges[b], edges[b+1] - 1
        i1 < i0 && continue
        seg = view(y, i0:i1)
        imin, imax = argmin(seg) + i0 - 1, argmax(seg) + i0 - 1
        for i in (min(imin, imax), max(imin, imax))
            push!(tt, t[i]); push!(yy, y[i])
        end
    end
    (tt, yy)
end

function shade_settle!(ax, r::SimResult; label = true, ytext = nothing)
    t0 = r.t_ch[1]
    t0 < 0 || return
    vspan!(ax, t0, 0.0; color = (S.WASH, 0.9))
    if label && ytext !== nothing
        text!(ax, t0, ytext; text = "settle", align = (:left, :top), offset = (4, -2), fontsize = 10.5, color = S.MUTED)
    end
end

function mark_deploy_end!(ax, p::SimParams, ytext)
    vlines!(ax, [p.T_deploy]; color = S.BASELINE, linewidth = 1)
    text!(ax, p.T_deploy, 0.0; text = "tip stops", align = (:left, :bottom), offset = (4, 3), fontsize = 10.5, color = S.MUTED)
end

segment_mid(lay::Layout) = 0.5 .* (lay.arc[1:end-1] .+ lay.arc[2:end])

function hinge_guides!(ax, lay::Layout)
    hlines!(ax, lay.arc[lay.hinge_nodes]; color = (S.INK, 0.12), linewidth = 0.75)
end

peak_time(r::SimResult) = r.peaks.F_tip >= r.peaks.F_root ? r.peaks.t_tip : r.peaks.t_root

# ------------------------------------------------------ 1. attachment loads --

function plot_attachment_forces(r::SimResult, path)
    p = r.params
    t = r.t_ch
    Ftip, Froot = r.ch["F_tip"], r.ch["F_root"]
    w = p.width
    pk = r.peaks
    tp = peak_time(r)
    fig = Figure(size = (1280, 760))
    header!(fig[1, 1:2], "Attachment loads during deployment",
        "Magnitude of the force each attachment applies to the structure", run_context(p))

    left = GridLayout(fig[2, 1])
    axa = Axis(left[1, 1]; ylabel = "Tip anchor\nposition [m]", height = 80,
        xticklabelsvisible = false, yticks = [0, 2, 4])
    shade_settle!(axa, r)
    lines!(axa, t, r.ch["anchor_x"]; color = S.INK2, linewidth = 1.5)
    xlims!(axa, t[1], t[end])

    ymax = 2.2 * max(pk.F_tip, pk.F_root, maximum(Ftip), maximum(Froot))
    axF = Axis(left[2, 1]; xlabel = "Time since deployment start [s]",
        ylabel = "Load per metre width [N/m]  (log scale)",
        yscale = Makie.pseudolog10, yticks = log_ticks(ymax))
    shade_settle!(axF, r; ytext = ymax)
    mark_deploy_end!(axF, p, ymax)
    tr, yr = envelope(t, Froot, 3000)
    tt, yt = envelope(t, Ftip, 3000)
    lines!(axF, tr, yr; color = S.ROOT, linewidth = 1.1, label = "Root attachment")
    lines!(axF, tt, yt; color = S.TIP, linewidth = 1.1, label = "Tip attachment")
    peak_marker!(axF, pk.t_tip, pk.F_tip, S.TIP)
    peak_marker!(axF, pk.t_root, pk.F_root, S.ROOT)
    ylims!(axF, 0, ymax)
    xlims!(axF, t[1], t[end])
    linkxaxes!(axa, axF)
    axislegend(axF; position = :lt, orientation = :horizontal, margin = (60, 0, 0, 6))
    txt = @sprintf("Peak tip %s N/m  (%s N total)\nt = %s s", fmt(pk.F_tip), fmt(pk.F_tip * w), fmt(pk.t_tip; sig = 4))
    text!(axF, pk.t_tip, pk.F_tip; text = txt, align = (:right, :bottom), offset = (-10, 6), fontsize = 11.5, color = S.INK)
    rowgap!(left, 6)

    # Zoom on the peak event
    win = 0.06
    sel = findall(x -> abs(x - tp) <= win, t)
    axZ = Axis(fig[2, 2]; title = "Peak event", subtitle = "±60 ms around the largest load, 0.5 ms sampling",
        xlabel = "Time relative to peak [ms]", ylabel = "Load per metre width [N/m]")
    tz = (t[sel] .- tp) .* 1e3
    lines!(axZ, tz, Froot[sel]; color = S.ROOT, linewidth = 1.6)
    lines!(axZ, tz, Ftip[sel]; color = S.TIP, linewidth = 1.6)
    ipk = argmax(max.(Ftip[sel], Froot[sel]))
    ytop = max(Ftip[sel][ipk], Froot[sel][ipk])
    peak_marker!(axZ, tz[ipk], ytop, Ftip[sel][ipk] >= Froot[sel][ipk] ? S.TIP : S.ROOT)
    text!(axZ, tz[ipk], ytop; text = fmt(ytop) * " N/m", align = (:left, :bottom), offset = (8, 2), fontsize = 11.5, color = S.INK)
    ylims!(axZ, 0, 1.18 * ytop)
    xlims!(axZ, -win * 1e3, win * 1e3)
    preload = p.EA * p.preload_strain
    hlines!(axZ, [preload]; color = S.MUTED, linewidth = 1, linestyle = :dash)
    text!(axZ, -win * 1e3, preload; text = "preload EA·ε = " * fmt(preload) * " N/m", align = (:left, :bottom),
        offset = (4, 2), fontsize = 10.5, color = S.MUTED)
    colsize!(fig.layout, 2, Relative(0.3))
    colgap!(fig.layout, 28)

    footer!(fig[3, 1:2], p; extra = @sprintf("Totals assume a %s m wide blanket.", fmt(w)))
    rowgap!(fig.layout, 1, 14)
    save(path, fig; px_per_unit = 2)
    path
end

# ------------------------------------------------------------- 2. tension --

function plot_tension(r::SimResult, path)
    p = r.params
    lay = p.layout
    t = r.t_ch
    pk = r.peaks
    fig = Figure(size = (1280, 820))
    header!(fig[1, 1:4], "Tension along the blanket",
        "Axial force per metre width (spring + damper; negative = compression)", run_context(p))

    Tmax = r.ch["tension_max"]
    ymax = 2.2 * max(maximum(Tmax), pk.tension)
    ax1 = Axis(fig[2, 1]; ylabel = "Max tension\n[N/m] (log)", yscale = Makie.pseudolog10,
        yticks = log_ticks(ymax), xticklabelsvisible = false, height = 170)
    shade_settle!(ax1, r; ytext = ymax)
    mark_deploy_end!(ax1, p, ymax)
    te, ye = envelope(t, Tmax, 3000)
    lines!(ax1, te, ye; color = S.BLUE, linewidth = 1.1)
    preload = p.EA * p.preload_strain
    hlines!(ax1, [preload]; color = S.MUTED, linewidth = 1, linestyle = :dash)
    text!(ax1, t[1], preload; text = "preload " * fmt(preload) * " N/m", align = (:left, :bottom), offset = (4, 1), fontsize = 10.5, color = S.MUTED)
    peak_marker!(ax1, pk.t_tension, pk.tension, S.BLUE)
    text!(ax1, pk.t_tension, pk.tension; text = @sprintf("%s N/m at %s m from root", fmt(pk.tension), fmt(pk.arc_tension)),
        align = (:right, :bottom), offset = (-10, 2), fontsize = 11.5, color = S.INK)
    ylims!(ax1, 0, ymax); xlims!(ax1, t[1], t[end])

    mids = segment_mid(lay)
    M = maximum(abs, r.tension)
    ax2 = Axis(fig[3, 1]; xlabel = "Time since deployment start [s]", ylabel = "Arc length from root [m]",
        yticks = arc_ticks(p))
    hm = heatmap!(ax2, r.t_f, mids, r.tension; colormap = S.DIVERGING, colorrange = (-M, M),
        colorscale = TENSION_SCALE)
    hinge_guides!(ax2, lay)
    linkxaxes!(ax1, ax2)
    xlims!(ax2, t[1], t[end])
    Colorbar(fig[3, 2], hm; label = "Tension [N/m]  (symmetric log; blue = compression)", ticks = sym_ticks(M))

    # High-rate window around the peak attachment load
    if !isempty(r.snap_t)
        tp = peak_time(r)
        ts = (r.snap_t .- tp) .* 1e3
        Ms = maximum(abs, r.snap_T)
        ax3 = Axis(fig[2:3, 3]; title = "Tension waves at the peak",
            subtitle = "Full-rate window (0.5 ms) around the largest attachment load",
            xlabel = "Time relative to peak [ms]", ylabel = "Arc length from root [m]",
            yticks = arc_ticks(p))
        hm3 = heatmap!(ax3, ts, mids, r.snap_T; colormap = S.DIVERGING, colorrange = (-Ms, Ms))
        hinge_guides!(ax3, lay)
        c = sqrt(p.EA / lay.areal_density)
        L = flat_length(p)
        tl = (L / c) * 1e3
        x0 = clamp(ts[1] + 0.12 * (ts[end] - ts[1]), ts[1], ts[end] - tl)
        lines!(ax3, [x0, x0 + tl], [L, 0.0]; color = (S.INK, 0.7), linewidth = 1, linestyle = :dash)
        text!(ax3, x0, L; text = @sprintf("axial wave transit\n√(EA/ρ) = %s m/s", fmt(c)),
            align = (:left, :top), offset = (10, -6), fontsize = 10.5, color = S.INK)
        xlims!(ax3, ts[1], ts[end]); ylims!(ax3, 0, L)
        Colorbar(fig[2:3, 4], hm3; label = "Tension [N/m]")
        colsize!(fig.layout, 3, Relative(0.3))
    end
    footer!(fig[4, 1:4], p)
    rowgap!(fig.layout, 1, 14)
    save(path, fig; px_per_unit = 2)
    path
end

"""Symmetric log color scale, linear within +-0.1 N/m (shows slack-phase tension)."""
const TENSION_SCALE = Makie.Symlog10(0.1)

function arc_ticks(p::SimParams)
    L = flat_length(p)
    vals = collect(0.0:1.0:floor(L))
    labels = fmt.(vals)
    labels[1] = "0 root"
    if L - vals[end] < 0.35
        vals[end] = L
        labels[end] = fmt(L) * " tip"
    else
        push!(vals, L); push!(labels, fmt(L) * " tip")
    end
    (vals, labels)
end

function sym_ticks(M)
    vals = [0.0]
    for e in -1:5
        v = 10.0^e
        v > M && break
        push!(vals, v); pushfirst!(vals, -v)
    end
    (vals, fmt.(vals))
end

# --------------------------------------------------------- 3. fold angles --

function plot_fold_angles(r::SimResult, path)
    p = r.params
    lay = p.layout
    nh = length(lay.hinge_nodes)
    nh == 0 && return nothing
    colors = S.ordinal(nh)
    rest = p.theta_rest
    fig = Figure(size = (1280, 160 + 88 * nh))
    header!(fig[1, 1], "Fold hinge angles",
        "Fold angle per hinge: π = stowed, 0 = flat (dashed: rest angle = fold memory)", run_context(p))
    g = GridLayout(fig[2, 1])
    axes = Axis[]
    tf = r.t_f
    for k in 1:nh
        last = k == nh
        side = k == 1 ? "  (root side)" : last ? "  (tip side)" : ""
        ax = Axis(g[k, 1]; yticks = ([0, pi], ["0", "π"]), xticklabelsvisible = last,
            xlabel = last ? "Time since deployment start [s]" : "", ygridvisible = false,
            bottomspinevisible = last)
        push!(axes, ax)
        t0 = tf[1]
        t0 < 0 && vspan!(ax, t0, 0; color = (S.WASH, 0.9))
        hlines!(ax, [0.0]; color = S.BASELINE, linewidth = 0.75)
        hlines!(ax, [rest]; color = S.MUTED, linewidth = 0.9, linestyle = :dash)
        phi = r.fold[:, k]
        band!(ax, tf, zeros(length(tf)), phi; color = (colors[k], 0.14))
        lines!(ax, tf, phi; color = colors[k], linewidth = 1.6)
        vlines!(ax, [p.T_deploy]; color = S.BASELINE, linewidth = 1)
        text!(ax, tf[1], pi; text = "Hinge $k" * side, align = (:left, :top), offset = (6, 2), fontsize = 11.5, color = S.INK)
        xlims!(ax, tf[1], tf[end])
        lo = min(-0.1pi, minimum(phi) - 0.05)
        hi = max(1.08pi, maximum(phi) + 0.05)
        ylims!(ax, lo, hi)
    end
    text!(axes[1], p.T_deploy, pi; text = "tip stops", align = (:right, :top), offset = (-4, 2), fontsize = 10.5, color = S.MUTED)
    text!(axes[end], tf[end], rest; text = @sprintf("rest %.2gπ", rest / pi), align = (:right, :bottom), offset = (-4, 1), fontsize = 10.5, color = S.MUTED)
    linkxaxes!(axes...)
    rowgap!(g, 4)
    footer!(fig[3, 1], p)
    rowgap!(fig.layout, 1, 14)
    save(path, fig; px_per_unit = 2)
    path
end

# ------------------------------------------------------ 4. cell curvature --

function plot_curvature(r::SimResult, path)
    p = r.params
    lay = p.layout
    pk = r.peaks
    t = r.t_ch
    ratio = r.ch["curvature_ratio_max"]
    fig = Figure(size = (1280, 800))
    header!(fig[1, 1:3], "Cell damage check: panel curvature",
        @sprintf("Curvature inside panels relative to the allowable 1/R_min (R_min = %s m); above 1 exceeds the cell bend limit", fmt(p.min_cell_bend_radius)), run_context(p))

    ymax = max(1.3, 1.12 * max(maximum(ratio), pk.curvature_ratio))
    ax1 = Axis(fig[2, 1]; ylabel = "Max κ · R_min", xticklabelsvisible = false, height = 190)
    shade_settle!(ax1, r; ytext = ymax)
    te, ye = envelope(t, ratio, 3000)
    over = max.(ye, 1.0)
    band!(ax1, te, fill(1.0, length(te)), over; color = (S.CRITICAL, 0.12))
    lines!(ax1, te, ye; color = S.BLUE, linewidth = 1.2)
    hlines!(ax1, [1.0]; color = S.CRITICAL, linewidth = 1.5)
    text!(ax1, t[end], 1.0; text = "allowable", align = (:right, :top), offset = (-4, -3), fontsize = 11, color = S.INK2)
    peak_marker!(ax1, pk.t_curvature, pk.curvature_ratio, S.BLUE)
    exceeds = pk.curvature_ratio > 1
    msg = @sprintf("%s× allowable at t = %s s, %s m from root", fmt(pk.curvature_ratio), fmt(pk.t_curvature), fmt(pk.arc_curvature))
    text!(ax1, pk.t_curvature, pk.curvature_ratio; text = msg, align = (:left, :bottom), offset = (10, 2), fontsize = 11.5, color = S.INK)
    ylims!(ax1, 0, ymax); xlims!(ax1, t[1], t[end])

    C = abs.(r.curv)
    cmax = max(1.0, maximum(C))
    ax2 = Axis(fig[3, 1]; xlabel = "Time since deployment start [s]", ylabel = "Arc length from root [m]",
        yticks = 0:1:ceil(Int, flat_length(p)))
    hm = heatmap!(ax2, r.t_f, lay.arc, C; colormap = S.SEQUENTIAL, colorrange = (0, cmax))
    if cmax > 1
        contour!(ax2, r.t_f, lay.arc, C; levels = [1.0], color = S.CRITICAL, linewidth = 1)
    end
    hinge_guides!(ax2, lay)
    linkxaxes!(ax1, ax2)
    xlims!(ax2, t[1], t[end])
    Colorbar(fig[3, 2], hm; label = "|κ| · R_min  (red contour = allowable)")

    # Peak ratio by panel (all saved frames)
    peaks_by_panel = zeros(lay.n_panels)
    for i in 1:lay.N
        lay.kind[i] == PANEL_NODE || continue
        peaks_by_panel[lay.panel[i]] = max(peaks_by_panel[lay.panel[i]], maximum(view(C, :, i)))
    end
    ax3 = Axis(fig[2:3, 3]; title = "Peak by panel", subtitle = "Saved frames (every $(fmt(p.field_dt)) s)",
        xlabel = "Max |κ| · R_min", yticks = (1:lay.n_panels, ["Panel $k" for k in 1:lay.n_panels]),
        yreversed = true, xgridvisible = true, ygridvisible = false, bottomspinevisible = true)
    cols = [v > 1 ? S.CRITICAL : S.BLUE for v in peaks_by_panel]
    barplot!(ax3, 1:lay.n_panels, peaks_by_panel; direction = :x, color = cols, gap = 0.45, strokewidth = 0)
    vlines!(ax3, [1.0]; color = S.INK2, linewidth = 1)
    xm = 1.7 * maximum(peaks_by_panel)
    for (k, v) in enumerate(peaks_by_panel)
        text!(ax3, v, k; text = fmt(v; sig = 2) * "×" * (v > 1 ? " exceeds" : ""), align = (:left, :center), offset = (6, 0), fontsize = 11, color = S.INK2)
    end
    xlims!(ax3, 0, max(xm, 1.4))
    text!(ax3, 1.0, lay.n_panels + 0.5; text = "allowable", align = (:left, :top), offset = (4, 0), fontsize = 10.5, color = S.MUTED)
    colsize!(fig.layout, 3, Relative(0.24))

    extra = exceeds ? "Red marks cell curvature above the placeholder allowable." : ""
    footer!(fig[4, 1:3], p; extra)
    rowgap!(fig.layout, 1, 14)
    save(path, fig; px_per_unit = 2)
    path
end

# ------------------------------------------------------------- 5. energy --

function plot_energy(r::SimResult, path)
    p = r.params
    t = r.t_ch
    ch = r.ch
    stored = ch["E_axial"] .+ ch["E_panel"] .+ ch["E_hinge"] .+ ch["E_att"]
    W = ch["W_tip"]
    Dkeys = ["D_axial", "D_panel", "D_hinge", "D_att", "D_drag", "D_settle"]
    Dnames = ["Axial damping", "Panel bending damping", "Hinge damping", "Attachment dampers", "Air drag", "Settle damping (numerical)"]
    D = sum(ch[k] for k in Dkeys)
    grav = maximum(abs, ch["PE_grav"] .- ch["PE_grav"][1]) > 1e-9

    fig = Figure(size = (1280, 780))
    header!(fig[1, 1:2], "Energy budget",
        "Work done by the tip attachment = Δ kinetic + Δ stored elastic + Δ gravitational + dissipated  (J per metre width)", run_context(p))

    ax = Axis(fig[2, 1]; ylabel = "Energy since run start [J/m]", xticklabelsvisible = false)
    shade_settle!(ax, r)
    series = [
        ("Work by tip attachment", W, S.INK),
        ("Dissipated", D, S.BLUE),
        ("Δ Stored elastic", stored .- stored[1], S.ORANGE),
        ("Kinetic", ch["KE"], S.AQUA),
    ]
    grav && push!(series, ("Δ Gravitational", ch["PE_grav"] .- ch["PE_grav"][1], S.YELLOW))
    for (name, y, col) in series
        te, ye = envelope(t, y, 2500)
        lines!(ax, te, ye; color = col, linewidth = name == "Work by tip attachment" ? 2.2 : 1.6, label = name)
    end
    vlines!(ax, [p.T_deploy]; color = S.BASELINE, linewidth = 1)
    lo, hi = extrema(vcat([y for (_, y, _) in series]...))
    span = hi - lo
    xlims!(ax, t[1], t[end])
    ylims!(ax, lo - 0.05span, hi + 0.08span)
    axislegend(ax; position = :lt, orientation = :horizontal, margin = (8, 0, 0, 6))

    axr = Axis(fig[3, 1]; xlabel = "Time since deployment start [s]", ylabel = "Residual\n[ppm of throughput]", height = 120)
    scale = maximum(abs, W) + maximum(D) + maximum(abs, stored .- stored[1])
    res = 1e6 .* ch["balance_residual"] ./ max(scale, eps())
    shade_settle!(axr, r)
    te, ye = envelope(t, res, 2000)
    lines!(axr, te, ye; color = S.INK2, linewidth = 1.2)
    rmax = max(maximum(abs, res), 1e-3)
    ylims!(axr, -1.6 * rmax, 1.6 * rmax)
    xlims!(axr, t[1], t[end])
    linkxaxes!(ax, axr)
    text!(axr, t[1], 1.6rmax; text = "Balance closes: max error " * sci(r.balance_error) * " of energy throughput (requirement < 1%)",
        align = (:left, :top), offset = (6, -4), fontsize = 11, color = S.INK)

    # Dissipation by mechanism (stacked)
    axd = Axis(fig[2:3, 2]; title = "Where the energy was dissipated", subtitle = "Cumulative, stacked",
        xlabel = "Time since deployment start [s]", ylabel = "J/m")
    base = zeros(length(t))
    used = [(k, n) for (k, n) in zip(Dkeys, Dnames) if maximum(abs, ch[k]) > 1e-9]
    stride = max(1, length(t) ÷ 2000)
    idx = 1:stride:length(t)
    for (j, (k, name)) in enumerate(used)
        top = base .+ ch[k]
        band!(axd, t[idx], base[idx], top[idx]; color = (S.SERIES[j], 0.85), label = name)
        base = top
    end
    xlims!(axd, t[1], t[end])
    ylims!(axd, 0, 1.05 * max(base[end], eps()))
    axislegend(axd; position = :lt, margin = (6, 0, 0, 6))
    colsize!(fig.layout, 2, Relative(0.32))
    footer!(fig[4, 1:2], p)
    rowgap!(fig.layout, 1, 14)
    save(path, fig; px_per_unit = 2)
    path
end

# ----------------------------------------------------------- 6. filmstrip --

"""Frame indices near the given times."""
frames_at(r, times) = unique([argmin(abs.(r.t_f .- τ)) for τ in times])

function blanket_segments(r::SimResult, k)
    N = size(r.pos, 2)
    pts = Vector{Point2f}(undef, 2(N - 1))
    for j in 1:N-1
        pts[2j-1] = Point2f(r.pos[k, j, 1], r.pos[k, j, 2])
        pts[2j] = Point2f(r.pos[k, j+1, 1], r.pos[k, j+1, 2])
    end
    pts
end

function seg_colors(r::SimResult, k)
    T = view(r.tension, k, :)
    repeat(collect(T); inner = 2)
end

function zlimits(r::SimResult)
    lo = minimum(view(r.pos, :, :, 2))
    hi = maximum(view(r.pos, :, :, 2))
    pad = 0.08 * max(hi - lo, 0.2)
    (lo - pad, hi + pad)
end

function xlimits(r::SimResult)
    L = flat_length(r.params)
    lo = min(-0.04L, minimum(view(r.pos, :, :, 1)) - 0.03L)
    hi = max(1.04L, maximum(view(r.pos, :, :, 1)) + 0.03L)
    (lo, hi)
end

"""Draw one blanket snapshot (true scale) into `ax`; returns the line plot."""
function draw_shape!(ax, r::SimResult, k, M; linewidth = 2.6, hinges = true, anchors = true)
    lay = r.params.layout
    hlines!(ax, [0.0]; color = S.GRID, linewidth = 0.75)
    ls = linesegments!(ax, blanket_segments(r, k); color = seg_colors(r, k), colormap = S.DIVERGING_LINE,
        colorrange = (-M, M), colorscale = TENSION_SCALE, linewidth)
    if hinges
        scatter!(ax, r.pos[k, lay.hinge_nodes, 1], r.pos[k, lay.hinge_nodes, 2]; color = S.SURFACE,
            strokecolor = S.INK2, strokewidth = 1, markersize = 5)
    end
    if anchors
        scatter!(ax, [r.pos[k, 1, 1]], [r.pos[k, 1, 2]]; marker = :rect, color = S.ROOT, markersize = 8)
        scatter!(ax, [r.pos[k, end, 1]], [r.pos[k, end, 2]]; marker = :rect, color = S.TIP, markersize = 8)
    end
    ls
end

function plot_filmstrip(r::SimResult, path; ncols = 3, nrows = 3)
    p = r.params
    tp = peak_time(r)
    n = ncols * nrows
    times = collect(range(0, p.T_deploy; length = n - 2))
    append!(times, [tp, p.T_deploy + p.T_hold])
    order = sortperm(times)
    ks = frames_at(r, times[order])
    while length(ks) < n                       # keep the grid full if frames coincide
        push!(ks, length(r.t_f))
    end
    zlo, zhi = zlimits(r)
    L = flat_length(p)
    M = maximum(abs, r.tension)
    cellw = 390
    xlo, xhi = xlimits(r)
    cellh = round(Int, cellw * (zhi - zlo) / (xhi - xlo))
    fig = Figure(size = (ncols * (cellw + 24) + 150, 150 + nrows * (cellh + 46)))
    header!(fig[1, 1:2], "Deployment sequence",
        "Blanket shape at true scale, colored by tension; open dots are fold hinges, squares are the root (orange) and tip (blue) attachment nodes",
        run_context(p))
    g = GridLayout(fig[2, 1])
    local ls
    for (j, k) in enumerate(ks[1:n])
        row, col = divrem(j - 1, ncols)
        ax = Axis(g[row+1, col+1]; aspect = DataAspect(), width = cellw, height = cellh,
            ygridvisible = false, yticksvisible = false, yticklabelsvisible = false,
            xticks = 0:1:floor(Int, L), xticklabelsize = 10.5)
        ls = draw_shape!(ax, r, k, M)
        τ = r.t_f[k]
        label = @sprintf("t = %s s", fmt(τ; sig = 4))
        note = abs(τ - tp) < 0.02 ? "peak attachment load" : τ > p.T_deploy + 1e-9 ? "tip stopped" :
            @sprintf("%.0f%% of tip travel", 100 * progress(p.profile, τ)[1])
        ax.title = label
        ax.subtitle = note
        ax.titlesize = 12.5
        xlims!(ax, xlimits(r)...); ylims!(ax, zlo, zhi)
    end
    colgap!(g, 24); rowgap!(g, 14)
    Colorbar(fig[2, 2], ls; label = "Tension [N/m]  (symmetric log; blue = compression)", ticks = sym_ticks(M),
        height = Relative(0.7))
    footer!(fig[3, 1:2], p)
    rowgap!(fig.layout, 1, 14)
    resize_to_layout!(fig)
    save(path, fig; px_per_unit = 2)
    path
end

# ------------------------------------------------------------ 7. summary --

function stat_tile!(pos, label, value, detail; status = nothing)
    Box(pos; color = S.PLANE, strokecolor = S.GRID, strokewidth = 1, cornerradius = 8)
    g = GridLayout(pos; alignmode = Outside(14))
    Label(g[1, 1], label; fontsize = 11.5, color = S.INK2, halign = :left, tellwidth = false)
    Label(g[2, 1], value; fontsize = 26, font = :medium, color = S.INK, halign = :left, tellwidth = false)
    if status === nothing
        Label(g[3, 1], detail; fontsize = 11, color = S.INK2, halign = :left, tellwidth = false)
    else
        gg = GridLayout(g[3, 1]; tellwidth = false, halign = :left)
        ax = Axis(gg[1, 1]; width = 10, height = 10, backgroundcolor = :transparent)
        hidedecorations!(ax); hidespines!(ax)
        marker, color = status == :bad ? (:utriangle, S.CRITICAL) : (:circle, S.GOOD)
        scatter!(ax, [0.0], [0.0]; marker, color, markersize = status == :bad ? 10 : 8)
        Label(gg[1, 2], detail; fontsize = 11, color = S.INK2, halign = :left)
        colgap!(gg, 5)
    end
    rowgap!(g, 2)
end

function plot_summary(r::SimResult, path)
    p = r.params
    pk = r.peaks
    w = p.width
    t = r.t_ch
    fig = Figure(size = (1280, 900))
    header!(fig[1, 1], "Flarewing blanket deployment: run summary", run_context(p))

    tiles = GridLayout(fig[2, 1])
    preload = p.EA * p.preload_strain
    stat_tile!(tiles[1, 1], "Peak tip load", fmt(pk.F_tip) * " N/m", @sprintf("%s N total · t = %s s", fmt(pk.F_tip * w), fmt(pk.t_tip; sig = 4)))
    stat_tile!(tiles[1, 2], "Peak root load", fmt(pk.F_root) * " N/m", @sprintf("%s N total · t = %s s", fmt(pk.F_root * w), fmt(pk.t_root; sig = 4)))
    stat_tile!(tiles[1, 3], "Peak blanket tension", fmt(pk.tension) * " N/m", @sprintf("%s× the %s N/m preload", fmt(pk.tension / preload; sig = 2), fmt(preload)))
    exceed = pk.curvature_ratio > 1
    stat_tile!(tiles[1, 4], "Max cell curvature", fmt(pk.curvature_ratio; sig = 2) * "× allowable",
        exceed ? @sprintf("Exceeds R_min = %s m", fmt(p.min_cell_bend_radius)) : "Within allowable"; status = exceed ? :bad : :good)
    ok = r.balance_error < 0.01
    stat_tile!(tiles[1, 5], "Energy balance error", sci(r.balance_error),
        ok ? "Closes (requirement < 1%)" : "Does not close"; status = ok ? :good : :bad)
    colgap!(tiles, 12)

    body = GridLayout(fig[3, 1])
    ymax = 2.2 * max(pk.F_tip, pk.F_root)
    axF = Axis(body[1, 1]; title = "Attachment loads", subtitle = "Per metre width, log scale",
        xlabel = "Time since deployment start [s]", ylabel = "N/m", yscale = Makie.pseudolog10, yticks = log_ticks(ymax))
    shade_settle!(axF, r; ytext = ymax)
    vlines!(axF, [p.T_deploy]; color = S.BASELINE, linewidth = 1)
    tr, yr = envelope(t, r.ch["F_root"], 1500)
    tt, yt = envelope(t, r.ch["F_tip"], 1500)
    lines!(axF, tr, yr; color = S.ROOT, linewidth = 1, label = "Root")
    lines!(axF, tt, yt; color = S.TIP, linewidth = 1, label = "Tip")
    peak_marker!(axF, pk.t_tip, pk.F_tip, S.TIP)
    ylims!(axF, 0, ymax); xlims!(axF, t[1], t[end])
    axislegend(axF; position = :lt, orientation = :horizontal, margin = (50, 0, 0, 4))

    lay = p.layout
    nh = length(lay.hinge_nodes)
    axH = Axis(body[2, 1]; title = "Fold angles", subtitle = "Each row is a hinge, root at top; light = flat, dark = stowed",
        xlabel = "Time since deployment start [s]", yticks = (1:nh, ["H$k" for k in 1:nh]), yreversed = true,
        ygridvisible = false)
    if nh > 0
        hmH = heatmap!(axH, r.t_f, 1:nh, clamp.(r.fold, 0, pi); colormap = S.SEQUENTIAL, colorrange = (0, pi))
        Colorbar(body[2, 2], hmH; ticks = ([0, pi / 2, pi], ["0", "π/2", "π"]), label = "fold angle")
    end
    xlims!(axH, t[1], t[end])
    linkxaxes!(axF, axH)

    # Mini filmstrip (true scale)
    tp = peak_time(r)
    times = sort([0.0, 0.2p.T_deploy, 0.45p.T_deploy, 0.7p.T_deploy, tp, p.T_deploy + p.T_hold])
    ks = frames_at(r, times)
    strip = GridLayout(body[1:2, 3])
    zlo, zhi = zlimits(r)
    M = maximum(abs, r.tension)
    L = flat_length(p)
    Label(strip[0, 1:2], "Blanket shape"; font = :medium, fontsize = 13.5, halign = :left, tellwidth = false)
    Label(strip[1, 1:2], "True scale, colored by tension (gray = slack, red = taut)"; fontsize = 11.5, color = S.INK2, halign = :left, tellwidth = false)
    for (j, k) in enumerate(ks)
        row, col = divrem(j - 1, 2)
        ax = Axis(strip[row+2, col+1]; aspect = DataAspect(), ygridvisible = false, yticksvisible = false,
            yticklabelsvisible = false, xticks = [0, L], xticklabelsize = 10, title = @sprintf("t = %s s", fmt(r.t_f[k]; sig = 3)),
            titlesize = 11.5, titlefont = :regular, titlecolor = S.INK2)
        draw_shape!(ax, r, k, M; linewidth = 2, hinges = false, anchors = false)
        xlims!(ax, xlimits(r)...); ylims!(ax, zlo, zhi)
    end
    rowgap!(strip, 6); colgap!(strip, 14)
    colsize!(body, 3, Relative(0.36))
    colgap!(body, 24)
    rowgap!(body, 16)
    footer!(fig[4, 1], p)
    rowgap!(fig.layout, 1, 14)
    rowgap!(fig.layout, 2, 18)
    save(path, fig; px_per_unit = 2)
    path
end

# ----------------------------------------------------------- 8. animation --

"""
    animate_run(r, dir; fps, seconds, gif) -> paths

MP4 (and optionally GIF) of the blanket shape colored by tension, with the
attachment-load trace revealed in sync. The deployment is time-compressed;
the second around the peak load plays in slow motion.
"""
function animate_run(r::SimResult, dir; fps = 30, seconds = 16.0, gif = true)
    p = r.params
    lay = p.layout
    tf = r.t_f
    tp = peak_time(r)
    t_first = max(tf[1], -0.5)
    nframes = round(Int, fps * seconds)
    # 25% of frames: +-0.5 s around the peak at the saved field rate; rest uniform.
    slow_lo, slow_hi = max(t_first, tp - 0.5), min(tf[end], tp + 0.5)
    slow_frames = findall(τ -> slow_lo <= τ <= slow_hi, tf)
    n_rest = max(nframes - length(slow_frames), 10)
    uniform = frames_at(r, range(t_first, tf[end]; length = n_rest))
    ks = sort(unique(vcat(filter(k -> !(slow_lo <= tf[k] <= slow_hi), uniform), slow_frames)))
    slowmo_factor = (slow_hi - slow_lo) / (length(slow_frames) / fps)

    M = maximum(abs, r.tension)
    zlo, zhi = zlimits(r)
    L = flat_length(p)

    fig = Figure(size = (1280, 720), figure_padding = (28, 28, 18, 18))
    top = GridLayout(fig[1, 1])
    Label(top[1, 1], "Flarewing blanket deployment"; fontsize = 22, font = :medium, halign = :left, tellwidth = false)
    Label(top[2, 1], run_context(p); fontsize = 12.5, color = S.INK2, halign = :left, tellwidth = false)
    time_txt = Observable("t = 0.00 s")
    phase_txt = Observable("")
    Label(top[1, 2], time_txt; fontsize = 22, font = :medium, halign = :right, tellwidth = true)
    Label(top[2, 2], phase_txt; fontsize = 12.5, color = S.INK2, halign = :right, tellwidth = true)
    rowgap!(top, 2)

    axg = Axis(fig[2, 1]; aspect = DataAspect(), xlabel = "x: deployment direction [m]", ylabel = "z [m]",
        ygridvisible = false)
    hlines!(axg, [0.0]; color = S.GRID, linewidth = 0.75)
    k0 = ks[1]
    segs = Observable(blanket_segments(r, k0))
    cols = Observable(seg_colors(r, k0))
    hinge_pts = Observable([Point2f(r.pos[k0, i, 1], r.pos[k0, i, 2]) for i in lay.hinge_nodes])
    tip_node = Observable(Point2f(r.pos[k0, end, 1], r.pos[k0, end, 2]))
    anchor = Observable(Point2f(anchor_x_at(r, tf[k0]), 0.0))
    link = @lift [$anchor, $tip_node]
    ls = linesegments!(axg, segs; color = cols, colormap = S.DIVERGING_LINE, colorrange = (-M, M),
        colorscale = TENSION_SCALE, linewidth = 3.4)
    scatter!(axg, hinge_pts; color = S.SURFACE, strokecolor = S.INK2, strokewidth = 1.2, markersize = 7)
    lines!(axg, link; color = S.TIP, linewidth = 1)
    scatter!(axg, [Point2f(r.pos[1, 1, 1], r.pos[1, 1, 2])]; marker = :rect, color = S.ROOT, markersize = 12)
    scatter!(axg, anchor; marker = :rect, color = S.TIP, markersize = 12)
    text!(axg, r.pos[1, 1, 1], r.pos[1, 1, 2]; text = "root", align = (:right, :bottom), offset = (-8, 4), fontsize = 11.5, color = S.INK2)
    text!(axg, anchor; text = "tip anchor", align = (:center, :bottom), offset = (0, 10), fontsize = 11.5, color = S.INK2)
    xlo, xhi = xlimits(r)
    xlims!(axg, min(xlo, -0.08L), max(xhi, 1.08L)); ylims!(axg, zlo, zhi)
    Colorbar(fig[2, 2], ls; label = "Tension [N/m]", ticks = sym_ticks(M), height = Relative(0.8))

    t = r.t_ch
    ymax = 2.2 * max(r.peaks.F_tip, r.peaks.F_root)
    axf = Axis(fig[3, 1]; height = 150, xlabel = "Time since deployment start [s]", ylabel = "Load [N/m] (log)",
        yscale = Makie.pseudolog10, yticks = log_ticks(ymax))
    sel = findall(τ -> τ >= t_first, t)
    tr, yr = envelope(t[sel], r.ch["F_root"][sel], 2500)
    tt, yt = envelope(t[sel], r.ch["F_tip"][sel], 2500)
    lines!(axf, tr, yr; color = (S.ROOT, 0.18), linewidth = 1)
    lines!(axf, tt, yt; color = (S.TIP, 0.18), linewidth = 1)
    nshow_r = Observable(1); nshow_t = Observable(1)
    lines!(axf, @lift(tr[1:$nshow_r]), @lift(yr[1:$nshow_r]); color = S.ROOT, linewidth = 1.2, label = "Root")
    lines!(axf, @lift(tt[1:$nshow_t]), @lift(yt[1:$nshow_t]); color = S.TIP, linewidth = 1.2, label = "Tip")
    cursor = Observable(t_first)
    vlines!(axf, @lift([$cursor]); color = S.INK, linewidth = 1)
    xlims!(axf, t_first, t[end]); ylims!(axf, 0, ymax)
    axislegend(axf; position = :lt, orientation = :horizontal)
    Label(fig[4, 1:2], DISCLAIMER; fontsize = 11, color = S.INK2, halign = :left, tellwidth = false)
    rowgap!(fig.layout, 8)

    mp4 = joinpath(dir, "deployment.mp4")
    record(fig, mp4, ks; framerate = fps) do k
        τ = tf[k]
        segs[] = blanket_segments(r, k)
        cols[] = seg_colors(r, k)
        hinge_pts[] = [Point2f(r.pos[k, i, 1], r.pos[k, i, 2]) for i in lay.hinge_nodes]
        tip_node[] = Point2f(r.pos[k, end, 1], r.pos[k, end, 2])
        anchor[] = Point2f(anchor_x_at(r, τ), anchor_z_at(r, τ))
        time_txt[] = @sprintf("t = %.2f s", τ)
        phase = τ < 0 ? "settling (tip held)" : τ <= p.T_deploy ? @sprintf("deploying: %.1f%% of tip travel", 100 * progress(p.profile, τ)[1]) : "holding at full deployment"
        if slow_lo <= τ <= slow_hi
            phase = @sprintf("slow motion %.2g×  ·  peak load event  ·  ", slowmo_factor) * phase
        end
        phase_txt[] = phase
        nshow_r[] = max(1, searchsortedlast(tr, τ))
        nshow_t[] = max(1, searchsortedlast(tt, τ))
        cursor[] = τ
    end
    paths = [mp4]
    if gif
        g = joinpath(dir, "deployment.gif")
        ff = Makie.FFMPEG_jll.ffmpeg()
        filt = "fps=15,scale=960:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=96:stats_mode=diff[p];[s1][p]paletteuse=dither=bayer:bayer_scale=4"
        run(pipeline(`$ff -y -loglevel error -i $mp4 -vf $filt $g`))
        push!(paths, g)
    end
    paths
end

function anchor_x_at(r::SimResult, τ)
    i = clamp(searchsortedfirst(r.t_ch, τ), 1, length(r.t_ch))
    r.ch["anchor_x"][i]
end

function anchor_z_at(r::SimResult, τ)
    i = clamp(searchsortedfirst(r.t_ch, τ), 1, length(r.t_ch))
    r.ch["anchor_z"][i]
end

# -------------------------------------------------------------- driver --

"""
    plot_run(r, dir; animation=true) -> Vector of written paths

Write every figure for one run into `dir`.
"""
function plot_run(r::SimResult, dir; animation = true, fps = 30, seconds = 16.0, gif = true)
    mkpath(dir)
    out = String[]
    with_theme(Style.theme()) do
        push!(out, plot_summary(r, joinpath(dir, "summary.png")))
        push!(out, plot_attachment_forces(r, joinpath(dir, "attachment_forces.png")))
        push!(out, plot_tension(r, joinpath(dir, "tension.png")))
        fa = plot_fold_angles(r, joinpath(dir, "fold_angles.png"))
        fa === nothing || push!(out, fa)
        push!(out, plot_curvature(r, joinpath(dir, "panel_curvature.png")))
        push!(out, plot_energy(r, joinpath(dir, "energy.png")))
        push!(out, plot_filmstrip(r, joinpath(dir, "deployment_filmstrip.png")))
        if animation
            append!(out, animate_run(r, dir; fps, seconds, gif))
        end
    end
    out
end

# -------------------------------------------------------------- sweeps --

const PROFILE_ORDER = ["smooth", "constant_then_stop", "trapezoid"]
profile_color(prof) = S.SERIES[something(findfirst(==(prof), PROFILE_ORDER), 4)]
profile_label(prof) = Dict("smooth" => "Smooth S-curve", "constant_then_stop" => "Constant speed, abrupt stop",
    "trapezoid" => "Trapezoid")[prof]

"""
    plot_sweep_metric(agg, spec, metric, ylabel, path; title)

Small multiples, one panel per swept parameter, shared y-axis. Lines are the
jitter-ensemble mean; the shaded band spans the ensemble min-max.
"""
function plot_sweep_metric(agg, spec::SweepSpec, metric, ylabel, path; title, description, threshold = nothing)
    groups = spec.groups
    ncols = 4
    nrows = cld(length(groups) + 1, ncols)
    fig = Figure(size = (1400, 170 + 300 * nrows))
    ctx = @sprintf("One parameter varied at a time around %s  ·  line = mean of %d runs with T_deploy × {%s}, band = min–max (flapping-phase scatter)",
        spec.base["run_name"], length(spec.jitter), join(fmt.(spec.jitter), ", "))
    header!(fig[1, 1], title, description, ctx)
    grid = GridLayout(fig[2, 1])
    ymax = 1.12 * maximum(a[metric*"_max"] for a in agg)
    axes = Axis[]
    for (gi, g) in enumerate(groups)
        row, col = divrem(gi - 1, ncols)
        pts_all = filter(a -> a["group"] == g["name"], agg)
        categorical = haskey(g, "cases")
        logx = get(g, "scale", "linear") == "log"
        base = baseline_value(g, spec)
        sub = categorical ? "" : "key: " * g["key"]
        ax = Axis(grid[row+1, col+1]; title = g["name"], subtitle = sub,
            ylabel = col == 0 ? ylabel : "", yticklabelsvisible = col == 0,
            xscale = logx ? log10 : identity)
        push!(axes, ax)
        if categorical
            ax.xticks = (1:length(g["cases"]), [string(c["label"]) for c in g["cases"]])
            ax.xticklabelrotation = pi / 7
        else
            vals = sort(unique([a["value"] for a in pts_all]))
            labels = [first(a["label"] for a in pts_all if a["value"] == v) for v in vals]
            ax.xticks = (vals, labels)
        end
        if isfinite(base)
            vlines!(ax, [base]; color = S.BASELINE, linewidth = 1)
            text!(ax, base, ymax; text = "baseline", align = (:left, :top), offset = (3, -2), fontsize = 10, color = S.MUTED)
        end
        threshold === nothing || hlines!(ax, [threshold]; color = S.CRITICAL, linewidth = 1.2)
        profs = [p for p in PROFILE_ORDER if any(a -> a["profile"] == p, pts_all)]
        for prof in profs
            pts = sort(filter(a -> a["profile"] == prof, pts_all); by = a -> a["value"])
            x = [a["value"] for a in pts]
            col_ = profile_color(prof)
            band!(ax, x, [a[metric*"_min"] for a in pts], [a[metric*"_max"] for a in pts]; color = (col_, 0.13))
            lines!(ax, x, [a[metric*"_mean"] for a in pts]; color = col_, linewidth = 2)
            scatter!(ax, x, [a[metric*"_mean"] for a in pts]; color = col_, markersize = 8, strokecolor = S.SURFACE, strokewidth = 2)
        end
        ylims!(ax, 0, ymax)
        if categorical
            xlims!(ax, 0.5, length(g["cases"]) + 0.5)
        elseif logx
            vals = [a["value"] for a in pts_all]
            xlims!(ax, minimum(vals) / 1.35, maximum(vals) * 1.35)
        else
            vals = [a["value"] for a in pts_all]
            pad = 0.08 * (maximum(vals) - minimum(vals))
            xlims!(ax, minimum(vals) - pad, maximum(vals) + pad)
        end
    end
    linkyaxes!(axes...)
    # legend cell
    lrow, lcol = divrem(length(groups), ncols)
    elems = [[PolyElement(color = (profile_color(p), 0.13)), LineElement(color = profile_color(p), linewidth = 2),
        MarkerElement(marker = :circle, color = profile_color(p), markersize = 8, strokecolor = S.SURFACE, strokewidth = 2)]
        for p in PROFILE_ORDER if any(a -> a["profile"] == p, agg)]
    names = [profile_label(p) for p in PROFILE_ORDER if any(a -> a["profile"] == p, agg)]
    leg = GridLayout(grid[lrow+1, lcol+1]; tellheight = false, tellwidth = false, valign = :center, halign = :left)
    Label(leg[1, 1], "Deployment profile"; font = :medium, fontsize = 13, halign = :left, tellwidth = false)
    Legend(leg[2, 1], elems, names; halign = :left, orientation = :vertical, rowgap = 6)
    Label(leg[3, 1], "Trapezoid runs only in the\ndeployment-time sweep."; fontsize = 11, color = S.INK2, halign = :left, justification = :left, tellwidth = false)
    for c in 1:ncols
        colsize!(grid, c, Relative(1 / ncols))
    end
    colgap!(grid, 22); rowgap!(grid, 26)
    p0 = SimParams(spec.base)
    footer!(fig[3, 1], p0)
    rowgap!(fig.layout, 1, 16)
    save(path, fig; px_per_unit = 2)
    path
end

"""
    plot_sensitivity(sens, agg, spec, path)

Ranking of how much each parameter moves the peak attachment load across its
swept range, per profile, on a log fold-change axis around the baseline.
"""
function plot_sensitivity(agg, spec::SweepSpec, path; metric = "peak_attachment_force_N_per_m")
    profs = [p for p in PROFILE_ORDER if p in spec.profiles]
    groups = spec.groups
    # order by the largest fold change over profiles
    span(g, prof) = begin
        pts = filter(a -> a["group"] == g["name"] && a["profile"] == prof, agg)
        isempty(pts) ? 1.0 : maximum(a[metric*"_mean"] for a in pts) / minimum(a[metric*"_mean"] for a in pts)
    end
    order = sortperm([maximum(span(g, p) for p in profs) for g in groups]; rev = true)
    G = groups[order]
    fig = Figure(size = (1400, 200 + 62 * length(G)))
    header!(fig[1, 1:length(profs)], "What moves the peak attachment load?",
        "Range of the peak load across each parameter's swept range, relative to the baseline run (log scale). Longest bar = measure first.",
        @sprintf("Gray band: flapping-phase noise floor (median min–max spread of the %d-run jitter ensembles). Effects inside it are not resolvable.", length(spec.jitter)))
    allratios = Float64[]
    for prof in profs, g in G
        pts = filter(a -> a["group"] == g["name"] && a["profile"] == prof, agg)
        b = baseline_mean(agg, g, spec, prof; metric)
        for a in pts
            push!(allratios, a[metric*"_mean"] / b)
        end
    end
    lo = min(0.9, minimum(allratios)) / 2.4
    hi = max(1.1, maximum(allratios)) * 2.2
    ticks = [t for t in (0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0) if lo <= t <= hi]
    for (j, prof) in enumerate(profs)
        spreads = Float64[]
        for a in agg
            a["profile"] == prof && push!(spreads, (a[metric*"_max"] - a[metric*"_min"]) / a[metric*"_mean"])
        end
        noise = isempty(spreads) ? 0.0 : sort(spreads)[cld(length(spreads), 2)]
        ax = Axis(fig[2, j]; title = profile_label(prof), subtitle = @sprintf("baseline peak %s N/m", fmt(baseline_mean(agg, G[1], spec, prof; metric))),
            xscale = log2, xticks = (ticks, [t >= 1 ? "×" * fmt(t) : "×" * fmt(t) for t in ticks]),
            xlabel = "Peak load relative to baseline", yticks = (1:length(G), [g["name"] for g in G]),
            yreversed = true, xgridvisible = true, ygridvisible = false, yticklabelsvisible = j == 1,
            yticklabelcolor = S.INK, yticklabelsize = 12.5)
        vspan!(ax, 1 / (1 + noise / 2), 1 + noise / 2; color = (S.INK, 0.07))
        vlines!(ax, [1.0]; color = S.INK2, linewidth = 1)
        for (k, g) in enumerate(G)
            pts = sort(filter(a -> a["group"] == g["name"] && a["profile"] == prof, agg); by = a -> a[metric*"_mean"])
            isempty(pts) && continue
            b = baseline_mean(agg, g, spec, prof; metric)
            rmin, rmax = pts[1][metric*"_mean"] / b, pts[end][metric*"_mean"] / b
            linesegments!(ax, [Point2f(rmin, k), Point2f(rmax, k)]; color = (profile_color(prof), 0.35), linewidth = 10)
            scatter!(ax, [rmin, rmax], [k, k]; color = profile_color(prof), markersize = 11, strokecolor = S.SURFACE, strokewidth = 2)
            text!(ax, rmin, k; text = pts[1]["label"], align = (:right, :center), offset = (-9, 0), fontsize = 10.5, color = S.INK2)
            text!(ax, rmax, k; text = pts[end]["label"] * @sprintf("   ×%s range", fmt(rmax / rmin; sig = 2)), align = (:left, :center), offset = (9, 0), fontsize = 10.5, color = S.INK2)
        end
        xlims!(ax, lo, hi)
        ylims!(ax, length(G) + 0.6, 0.4)
    end
    colgap!(fig.layout, 30)
    footer!(fig[3, 1:length(profs)], SimParams(spec.base); extra = "Labels at the bar ends give the parameter value that produced that extreme.")
    rowgap!(fig.layout, 1, 16)
    save(path, fig; px_per_unit = 2)
    path
end

"""Write all sweep figures into `dir`."""
function plot_sweep(agg, spec::SweepSpec, dir)
    mkpath(dir)
    with_theme(Style.theme()) do
        p0 = SimParams(spec.base)
        [plot_sensitivity(agg, spec, joinpath(dir, "sensitivity.png")),
            plot_sweep_metric(agg, spec, "peak_attachment_force_N_per_m", "Peak load [N/m]",
                joinpath(dir, "sweep_peak_attachment_force.png"); title = "Peak attachment load across the parameter sweeps",
                description = "Largest of the root and tip attachment load magnitudes, per metre of width"),
            plot_sweep_metric(agg, spec, "peak_tension_N_per_m", "Peak tension [N/m]",
                joinpath(dir, "sweep_peak_tension.png"); title = "Peak blanket tension across the parameter sweeps",
                description = "Largest axial force anywhere in the blanket, per metre of width"),
            plot_sweep_metric(agg, spec, "max_panel_curvature_ratio", "Max κ · R_min",
                joinpath(dir, "sweep_curvature.png"); title = "Cell curvature check across the parameter sweeps",
                description = @sprintf("Peak panel curvature relative to the allowable 1/R_min (R_min = %s m); red line = allowable", fmt(p0.min_cell_bend_radius)),
                threshold = 1.0)]
    end
end
