# Parameter sweeps: one-at-a-time variations around a base config, run on all
# threads, summarized to CSV and plotted, plus a sensitivity ranking.
#
# Peak loads in lightly damped runs depend on the phase of residual blanket
# flapping when the blanket goes taut, so every sweep point can be run as a
# small ensemble with slightly jittered deployment times. The spread of that
# ensemble is the "phase noise floor": parameter effects smaller than it are
# not resolvable from single runs.

"""One run in a sweep."""
struct SweepJob
    group::String
    label::String          # human-readable parameter value
    value::Float64         # numeric value used for plotting (index for categorical)
    profile::String
    jitter::Float64
    cfg::Dict{String,Any}
end

"""
    SweepSpec

Parsed sweep YAML: base config file, profiles, jitter factors, resolution rule
and the list of groups (`key` + `values`, or categorical `cases`).
"""
struct SweepSpec
    base::Dict{String,Any}
    name::String
    profiles::Vector{String}
    jitter::Vector{Float64}
    segment_length::Float64          # target segment length (m); keeps resolution fixed
    groups::Vector{Dict{String,Any}}
    out_dir::String
end

function load_sweep(path::AbstractString)
    raw = YAML.load_file(path; dicttype = Dict{String,Any})
    base = load_config(joinpath(dirname(path), raw["base"]))
    get(raw, "overrides", nothing) === nothing || merge_config!(base, raw["overrides"])
    p0 = SimParams(base)
    seg = p0.layout.panel_length / p0.layout.nodes_per_panel
    SweepSpec(base, get(raw, "name", splitext(basename(path))[1]),
        string.(get(raw, "profiles", ["smooth"])),
        Float64.(get(raw, "jitter", [1.0])),
        Float64(get(raw, "segment_length", seg)),
        raw["groups"], get(raw, "out_dir", joinpath("results", "sweep")))
end

"""Apply the fixed-resolution rule: nodes per panel from the target segment length."""
function apply_resolution!(cfg, seg)
    Lp = num(cfg, "blanket.panel_length")
    set_param!(cfg, "blanket.nodes_per_panel", max(4, round(Int, Lp / seg)))
    cfg
end

function sweep_jobs(spec::SweepSpec)
    jobs = SweepJob[]
    L_total = num(spec.base, "blanket.n_panels") * num(spec.base, "blanket.panel_length")
    T0 = num(spec.base, "deployment.T_deploy")
    for g in spec.groups
        gname = g["name"]
        profiles = string.(get(g, "profiles", spec.profiles))
        cases = Tuple{String,Float64,Dict{String,Any}}[]
        if haskey(g, "cases")                         # categorical group
            for (i, c) in enumerate(g["cases"])
                push!(cases, (string(c["label"]), Float64(i), c["set"]))
            end
        else
            key = g["key"]
            for v in g["values"]
                x = parse_number(v, key)
                sets = Dict{String,Any}(key => x)
                if get(g, "fixed_total_length", false) && key == "blanket.n_panels"
                    sets[key] = round(Int, x)
                    sets["blanket.panel_length"] = L_total / x
                end
                label = key == "hinge.theta_rest" ? @sprintf("%.3gπ", x / pi) : fmt_sweep(x)
                push!(cases, (label, x, sets))
            end
        end
        for (label, x, sets) in cases, prof in profiles, j in spec.jitter
            cfg = deepcopy(spec.base)
            for (k, v) in sets
                set_param!(cfg, k, v)
            end
            set_param!(cfg, "deployment.profile", prof)
            T = num(cfg, "deployment.T_deploy")
            set_param!(cfg, "deployment.T_deploy", T * j)
            apply_resolution!(cfg, spec.segment_length)
            cfg["run_name"] = @sprintf("%s | %s = %s | %s | jitter %.2f", spec.name, gname, label, prof, j)
            push!(jobs, SweepJob(gname, label, x, prof, j, cfg))
        end
    end
    jobs
end

fmt_sweep(x) = abs(x) >= 1e4 || (abs(x) < 1e-2 && x != 0) ? @sprintf("%.3g", x) : @sprintf("%g", round(x; sigdigits = 3))

const SWEEP_COLUMNS = ["group", "label", "value", "profile", "jitter",
    "peak_tip_force_N_per_m", "peak_root_force_N_per_m", "peak_attachment_force_N_per_m",
    "peak_tension_N_per_m", "max_panel_curvature_ratio", "time_of_peak_s",
    "KE_at_deploy_end_J_per_m", "final_tip_force_N_per_m", "energy_balance_error",
    "nodes_per_panel", "dt_s", "wall_time_s", "diverged"]

function run_job(job::SweepJob)
    p = SimParams(job.cfg)
    r = simulate(p; verbose = false, record_fields = false)
    s = summary_metrics(r)
    i_end = clamp(searchsortedfirst(r.t_ch, p.T_deploy), 1, length(r.t_ch))
    Dict{String,Any}(
        "group" => job.group, "label" => job.label, "value" => job.value, "profile" => job.profile,
        "jitter" => job.jitter,
        "peak_tip_force_N_per_m" => s["peak_tip_force_N_per_m"],
        "peak_root_force_N_per_m" => s["peak_root_force_N_per_m"],
        "peak_attachment_force_N_per_m" => max(s["peak_tip_force_N_per_m"], s["peak_root_force_N_per_m"]),
        "peak_tension_N_per_m" => s["peak_tension_N_per_m"],
        "max_panel_curvature_ratio" => s["max_panel_curvature_ratio"],
        "time_of_peak_s" => s["time_of_peak_attachment_force_s"],
        "KE_at_deploy_end_J_per_m" => r.ch["KE"][i_end],
        "final_tip_force_N_per_m" => s["final_tip_force_N_per_m"],
        "energy_balance_error" => s["energy_balance_error"],
        "nodes_per_panel" => p.layout.nodes_per_panel, "dt_s" => p.dt,
        "wall_time_s" => r.wall_time, "diverged" => r.diverged)
end

"""
    run_sweep(spec) -> Vector{Dict}

Run every job on all available threads (start Julia with `-t auto`). Jobs that
share a configuration (the baseline appears in every group) are run once.
"""
function run_sweep(spec::SweepSpec; io = stdout, cache_dir = joinpath(spec.out_dir, "runs"))
    jobs = sweep_jobs(spec)
    keyof(job) = YAML.write(merge(job.cfg, Dict("run_name" => "")))
    unique_keys = unique(keyof.(jobs))
    index = Dict(k => i for (i, k) in enumerate(unique_keys))
    first_job = Dict{String,SweepJob}()
    for job in jobs
        get!(first_job, keyof(job), job)
    end
    # Resumable: each finished run is cached as JSON keyed by a hash of its full config.
    mkpath(cache_dir)
    cache_file(k) = joinpath(cache_dir, string(hash(k); base = 16) * ".json")
    results = Vector{Dict{String,Any}}(undef, length(unique_keys))
    todo = Int[]
    for (i, k) in enumerate(unique_keys)
        f = cache_file(k)
        if isfile(f)
            results[i] = Dict{String,Any}(JSON.parsefile(f))
        else
            push!(todo, i)
        end
    end
    @printf(io, "[sweep %s] %d runs (%d unique, %d cached) on %d threads\n", spec.name, length(jobs),
        length(unique_keys), length(unique_keys) - length(todo), Threads.nthreads())
    flush(io)
    # Longest runs first so the thread pool stays busy at the end.
    cost(i) = (j = first_job[unique_keys[i]]; num(j.cfg, "deployment.T_deploy") / num(j.cfg, "blanket.panel_length") * num(j.cfg, "blanket.nodes_per_panel"))
    sort!(todo; by = cost, rev = true)
    done = Threads.Atomic{Int}(0)
    t0 = time()
    lk = ReentrantLock()
    # Shared work queue: every thread pulls the next job, so no thread idles
    # while another still holds a backlog (fixed chunking would).
    queue = Channel{Int}(length(todo))
    foreach(i -> put!(queue, i), todo)
    close(queue)
    workers = map(1:Threads.nthreads()) do _
        Threads.@spawn for i in queue
            res = run_job(first_job[unique_keys[i]])
            results[i] = res
            lock(lk) do
                n = Threads.atomic_add!(done, 1) + 1
                open(f -> JSON.print(f, res), cache_file(unique_keys[i]), "w")
                @printf(io, "  %3d/%d done  (%.0f s elapsed)\n", n, length(todo), time() - t0)
                flush(io)
            end
        end
    end
    foreach(wait, workers)
    rows = Dict{String,Any}[]
    for job in jobs
        row = copy(results[index[keyof(job)]])
        row["group"], row["label"], row["value"], row["profile"], row["jitter"] = job.group, job.label, job.value, job.profile, job.jitter
        push!(rows, row)
    end
    rows
end

"""Read `sweep_runs.csv` back into row dictionaries (numbers parsed)."""
function read_sweep_csv(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    rows = Dict{String,Any}[]
    for line in lines[2:end]
        cells = String[]
        buf = IOBuffer(); inq = false
        for ch in line
            if ch == '"'
                inq = !inq
            elseif ch == ',' && !inq
                push!(cells, String(take!(buf)))
            else
                write(buf, ch)
            end
        end
        push!(cells, String(take!(buf)))
        row = Dict{String,Any}()
        for (h, c) in zip(header, cells)
            x = tryparse(Float64, c)
            row[h] = h in ("group", "label", "profile") || x === nothing ? c : x
        end
        push!(rows, row)
    end
    rows
end

function write_csv(path, rows, columns)
    open(path, "w") do io
        println(io, join(columns, ","))
        for r in rows
            println(io, join([csv_cell(r[c]) for c in columns], ","))
        end
    end
    path
end

csv_cell(x::AbstractString) = occursin(',', x) || occursin('"', x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : x
csv_cell(x::Real) = string(x)
csv_cell(x) = string(x)

"""Aggregate jitter ensembles: mean/min/max of each metric per (group, label, profile)."""
function aggregate(rows)
    metrics = ["peak_attachment_force_N_per_m", "peak_tip_force_N_per_m", "peak_tension_N_per_m",
        "max_panel_curvature_ratio", "KE_at_deploy_end_J_per_m", "final_tip_force_N_per_m"]
    groups = Dict{Tuple{String,String,String},Vector{Dict{String,Any}}}()
    order = Tuple{String,String,String}[]
    for r in rows
        k = (r["group"], r["label"], r["profile"])
        haskey(groups, k) || push!(order, k)
        push!(get!(groups, k, Dict{String,Any}[]), r)
    end
    out = Dict{String,Any}[]
    for k in order
        rs = groups[k]
        a = Dict{String,Any}("group" => k[1], "label" => k[2], "profile" => k[3], "value" => rs[1]["value"], "n" => length(rs))
        for m in metrics
            v = [Float64(r[m]) for r in rs]
            a[m*"_mean"] = sum(v) / length(v)
            a[m*"_min"] = minimum(v)
            a[m*"_max"] = maximum(v)
        end
        push!(out, a)
    end
    out
end

"""
    sensitivity_table(agg, spec) -> Vector{Dict}

For each (group, profile): fold change of the mean peak attachment load across
the swept range (max/min), the baseline value, and the phase-noise floor
(median relative spread of the jitter ensembles in that group).
"""
function sensitivity_table(agg, spec::SweepSpec; metric = "peak_attachment_force_N_per_m")
    out = Dict{String,Any}[]
    keys_ = unique([(a["group"], a["profile"]) for a in agg])
    for (g, prof) in keys_
        pts = filter(a -> a["group"] == g && a["profile"] == prof, agg)
        means = [a[metric*"_mean"] for a in pts]
        spreads = [(a[metric*"_max"] - a[metric*"_min"]) / max(a[metric*"_mean"], eps()) for a in pts]
        imin, imax = argmin(means), argmax(means)
        push!(out, Dict{String,Any}(
            "group" => g, "profile" => prof,
            "fold_change" => means[imax] / max(means[imin], eps()),
            "min_mean" => means[imin], "min_at" => pts[imin]["label"],
            "max_mean" => means[imax], "max_at" => pts[imax]["label"],
            "noise_floor" => sort(spreads)[cld(length(spreads), 2)],
            "labels" => [a["label"] for a in pts], "means" => means,
        ))
    end
    sort!(out; by = r -> -r["fold_change"])
    out
end

"""Baseline value of a sweep group (numeric groups) or the index of the baseline case."""
function baseline_value(g::Dict, spec::SweepSpec)
    if haskey(g, "cases")
        i = findfirst(c -> get(c, "baseline", false) == true, g["cases"])
        return i === nothing ? NaN : Float64(i)
    end
    num(spec.base, g["key"])
end

"""Mean of `metric` at the baseline point of group `g` for `profile` (NaN if absent)."""
function baseline_mean(agg, g::Dict, spec::SweepSpec, profile; metric = "peak_attachment_force_N_per_m")
    b = baseline_value(g, spec)
    i = findfirst(a -> a["group"] == g["name"] && a["profile"] == profile && isapprox(a["value"], b; rtol = 1e-9), agg)
    i === nothing ? NaN : agg[i][metric*"_mean"]
end
