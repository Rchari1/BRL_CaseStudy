# Results files and summary metrics.

"""
    summary_metrics(r::SimResult) -> Dict{String,Any}

Headline numbers for one run. Forces are per metre of width unless the key
ends in `_total_N` (scaled by `blanket_width`). All are placeholder-parameter
estimates.
"""
function summary_metrics(r::SimResult)
    p = r.params
    pk = r.peaks
    w = p.width
    t_peak = pk.F_tip >= pk.F_root ? pk.t_tip : pk.t_root
    Dict{String,Any}(
        "run_name" => p.name,
        "profile" => profile_name(p.profile),
        "T_deploy_s" => p.T_deploy,
        "n_panels" => p.layout.n_panels,
        "peak_root_force_N_per_m" => pk.F_root,
        "peak_root_force_total_N" => pk.F_root * w,
        "time_of_peak_root_s" => pk.t_root,
        "peak_tip_force_N_per_m" => pk.F_tip,
        "peak_tip_force_total_N" => pk.F_tip * w,
        "time_of_peak_tip_s" => pk.t_tip,
        "time_of_peak_attachment_force_s" => t_peak,
        "peak_tension_N_per_m" => pk.tension,
        "time_of_peak_tension_s" => pk.t_tension,
        "arc_of_peak_tension_m" => pk.arc_tension,
        "min_axial_force_N_per_m" => pk.min_tension,
        "max_panel_curvature_ratio" => pk.curvature_ratio,
        "time_of_max_curvature_s" => pk.t_curvature,
        "arc_of_max_curvature_m" => pk.arc_curvature,
        "cell_bend_limit_exceeded" => pk.curvature_ratio > 1,
        "energy_balance_error" => r.balance_error,
        "final_tip_force_N_per_m" => r.ch["F_tip"][end],
        "blanket_width_m" => w,
        "dt_s" => p.dt,
        "axial_stiffness_scale" => p.axial_stiffness_scale,
        "steps" => r.steps,
        "wall_time_s" => r.wall_time,
        "diverged" => r.diverged,
    )
end

"""One-line run summary (placeholder-parameter estimate)."""
function summary_line(r::SimResult)
    s = summary_metrics(r)
    warn = r.params.axial_stiffness_scale != 1.0 ? @sprintf(" [EA scaled x%.3g: peaks under-predicted]", r.params.axial_stiffness_scale) : ""
    @sprintf("%s | peak root %.4g N/m (%.4g N total) | peak tip %.4g N/m (%.4g N total) | at t = %.4g s | peak tension %.4g N/m | max panel curvature ratio %.3g%s | energy balance error %.2e | PLACEHOLDER PARAMETERS, NOT VALIDATED%s",
        s["run_name"], s["peak_root_force_N_per_m"], s["peak_root_force_total_N"],
        s["peak_tip_force_N_per_m"], s["peak_tip_force_total_N"], s["time_of_peak_attachment_force_s"],
        s["peak_tension_N_per_m"], s["max_panel_curvature_ratio"],
        s["cell_bend_limit_exceeded"] ? " (EXCEEDS cell bend limit)" : "",
        s["energy_balance_error"], warn)
end

"""
    save_results(dir, r) -> path

Write `results.jld2` (HDF5-based; readable from Python with h5py) holding
channels, fields, the snap window, summary metrics and the full config.
"""
function save_results(dir::AbstractString, r::SimResult)
    mkpath(dir)
    path = joinpath(dir, "results.jld2")
    p = r.params
    lay = p.layout
    jldopen(path, "w") do f
        f["README"] = "BlanketSim results. Placeholder parameters, not validated against test data. Forces per metre of width (N/m) unless noted; energies J per metre of width."
        f["config_yaml"] = YAML.write(p.cfg)
        f["summary_json"] = JSON.json(summary_metrics(r))
        f["channels/t"] = r.t_ch
        for (k, v) in r.ch
            f["channels/$k"] = v
        end
        f["fields/t"] = r.t_f
        f["fields/x"] = r.pos[:, :, 1]
        f["fields/z"] = r.pos[:, :, 2]
        f["fields/tension"] = r.tension
        f["fields/fold_angle"] = r.fold
        f["fields/curvature_ratio"] = r.curv
        f["snap/t"] = r.snap_t
        f["snap/tension"] = r.snap_T
        f["snap/F_root"] = r.snap_Froot
        f["snap/F_tip"] = r.snap_Ftip
        f["layout/arc"] = lay.arc
        f["layout/kind"] = Int.(lay.kind)
        f["layout/panel"] = lay.panel
        f["layout/hinge_nodes"] = lay.hinge_nodes
        f["layout/mass"] = lay.mass
    end
    path
end

"""Write summary.txt (the one-line summary) and summary.json."""
function write_summary(dir::AbstractString, r::SimResult)
    mkpath(dir)
    open(joinpath(dir, "summary.txt"), "w") do io
        println(io, summary_line(r))
    end
    open(joinpath(dir, "summary.json"), "w") do io
        JSON.print(io, summary_metrics(r), 2)
    end
end
