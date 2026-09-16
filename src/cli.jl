# Command-line entry points (see bin/blanket_sim.jl):
#
#   julia --project -t auto bin/blanket_sim.jl run   configs/baseline.yaml [--out DIR] [--no-animation] [--set key=value ...]
#   julia --project -t auto bin/blanket_sim.jl sweep configs/sweep.yaml   [--out DIR]

const USAGE = """
usage:
  blanket_sim run   <config.yaml> [--out DIR] [--no-animation] [--set section.key=value ...]
  blanket_sim sweep <sweep.yaml>  [--out DIR] [--replot]   (--replot: redraw figures from sweep_runs.csv)
"""

function parse_flags(args)
    opts = Dict{String,Any}("set" => String[])
    pos = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out"
            opts["out"] = args[i+1]; i += 1
        elseif a == "--no-animation"
            opts["no-animation"] = true
        elseif a == "--replot"
            opts["replot"] = true
        elseif a == "--set"
            push!(opts["set"], args[i+1]); i += 1
        elseif startswith(a, "--")
            error("unknown option $a\n$USAGE")
        else
            push!(pos, a)
        end
        i += 1
    end
    pos, opts
end

"""Run one configuration and write results, summary, figures and animation."""
function cli_run(path::AbstractString; out = nothing, animation = true, sets = String[], io = stdout)
    cfg = load_config(path)
    for s in sets
        k, v = split(s, "="; limit = 2)
        set_param!(cfg, strip(k), YAML.load(v))
    end
    p = SimParams(cfg)
    dir = out === nothing ? joinpath(string(get_param(cfg, "output.dir")), p.name) : out
    mkpath(dir)
    r = simulate(p; io)
    r.diverged && @warn "Run diverged; outputs are incomplete."
    files = [save_results(dir, r)]
    write_summary(dir, r)
    push!(files, joinpath(dir, "summary.txt"), joinpath(dir, "summary.json"))
    println(io, "  rendering figures...")
    animate = animation && get_param(cfg, "output.animation") == true
    append!(files, plot_run(r, dir; animation = animate,
        fps = Int(num(cfg, "output.animation_fps")), seconds = num(cfg, "output.animation_seconds"),
        gif = get_param(cfg, "output.gif") == true))
    println(io, "\n", summary_line(r))
    println(io, "\nwrote:")
    foreach(f -> println(io, "  ", f), files)
    r
end

"""Run a sweep and write the CSVs, figures and a sensitivity table."""
function cli_sweep(path::AbstractString; out = nothing, replot = false, io = stdout)
    spec = load_sweep(path)
    dir = out === nothing ? spec.out_dir : out
    mkpath(dir)
    if replot
        rows = read_sweep_csv(joinpath(dir, "sweep_runs.csv"))
        println(io, "[sweep $(spec.name)] re-plotting $(length(rows)) runs from sweep_runs.csv")
    else
        rows = run_sweep(spec; io)
        write_csv(joinpath(dir, "sweep_runs.csv"), rows, SWEEP_COLUMNS)
    end
    agg = aggregate(rows)
    aggcols = ["group", "label", "value", "profile", "n"]
    for m in ["peak_attachment_force_N_per_m", "peak_tip_force_N_per_m", "peak_tension_N_per_m",
        "max_panel_curvature_ratio", "KE_at_deploy_end_J_per_m", "final_tip_force_N_per_m"], s in ("mean", "min", "max")
        push!(aggcols, m * "_" * s)
    end
    write_csv(joinpath(dir, "sweep_summary.csv"), agg, aggcols)
    sens = sensitivity_table(agg, spec)
    write_sensitivity_table(joinpath(dir, "sensitivity_table.md"), sens, agg, spec)
    figs = plot_sweep(agg, spec, dir)
    println(io, "\nSensitivity of peak attachment load (fold change across swept range):")
    for s in sens
        @printf(io, "  %-40s %-20s ×%.2f  (%s N/m at %s  →  %s N/m at %s; noise floor ±%.0f%%)\n",
            s["group"], s["profile"], s["fold_change"], fmt(s["min_mean"]), s["min_at"], fmt(s["max_mean"]), s["max_at"], 50s["noise_floor"])
    end
    println(io, "\nwrote:\n  ", joinpath(dir, "sweep_runs.csv"), "\n  ", joinpath(dir, "sweep_summary.csv"), "\n  ",
        joinpath(dir, "sensitivity_table.md"))
    foreach(f -> println(io, "  ", f), figs)
    rows, agg, sens
end

function write_sensitivity_table(path, sens, agg, spec::SweepSpec)
    open(path, "w") do io
        println(io, "# Sensitivity table (auto-generated)\n")
        println(io, "Peak attachment load (max of root and tip, N per metre width). Mean of the jitter ensemble ",
            "(T_deploy × ", join(fmt.(spec.jitter), ", "), "). Placeholder parameters; not validated.\n")
        println(io, "| Rank | Parameter | Profile | Fold change | Lowest mean (at) | Highest mean (at) | Phase noise (± of mean) |")
        println(io, "|---:|---|---|---:|---|---|---:|")
        for (i, s) in enumerate(sens)
            @printf(io, "| %d | %s | %s | ×%.2f | %s N/m (%s) | %s N/m (%s) | ±%.0f%% |\n", i, s["group"], s["profile"],
                s["fold_change"], fmt(s["min_mean"]), s["min_at"], fmt(s["max_mean"]), s["max_at"], 50s["noise_floor"])
        end
    end
    path
end

"""Entry point used by bin/blanket_sim.jl. Returns a process exit code."""
function main(args = ARGS)
    isempty(args) && (print(USAGE); return 1)
    cmd = args[1]
    pos, opts = parse_flags(args[2:end])
    if cmd == "run" && length(pos) == 1
        cli_run(pos[1]; out = get(opts, "out", nothing), animation = !get(opts, "no-animation", false), sets = opts["set"])
    elseif cmd == "sweep" && length(pos) == 1
        cli_sweep(pos[1]; out = get(opts, "out", nothing), replot = get(opts, "replot", false))
    else
        print(USAGE)
        return 1
    end
    0
end
