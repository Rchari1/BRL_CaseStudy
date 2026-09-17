# Build the standalone interactive explorer: docs/index.html
#
#   julia --project=. tools/build_explorer.jl [case ...]
#
# Reads results/<case>/results.jld2 for each case (run the CLI first; the .jld2
# files are git-ignored) and results/sweep/sweep_summary.csv, compresses shape
# frames, tension, fold angles and load channels into base64 typed arrays, and
# inlines them into tools/explorer_template.html. The output is one HTML file
# that opens offline in any browser. Placeholder-parameter results.

using JLD2, JSON, Base64, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_CASES = ["baseline", "snap_worst_case", "ground_test", "hysteretic_hinge"]

b64f32(x) = base64encode(reinterpret(UInt8, Float32.(vec(x))))
b64i16(x) = base64encode(reinterpret(UInt8, Int16.(vec(x))))

"""Indices keeping the min and max of `y` in each of `nbins` bins (peaks survive decimation)."""
function envelope_indices(y, nbins)
    n = length(y)
    edges = round.(Int, range(1, n + 1; length = nbins + 1))
    idx = Int[]
    for b in 1:nbins
        i0, i1 = edges[b], edges[b+1] - 1
        i1 < i0 && continue
        seg = view(y, i0:i1)
        append!(idx, sort([argmin(seg) + i0 - 1, argmax(seg) + i0 - 1]))
    end
    idx
end

function export_case(name; frame_stride = 5, node_stride = 2)
    path = joinpath(ROOT, "results", name, "results.jld2")
    isfile(path) || error("missing $path: run `julia --project=. bin/blanket_sim.jl run configs/$name.yaml` first")
    jldopen(path, "r") do f
        t = f["channels/t"]
        ch(k) = f["channels/$k"]
        idx = sort(unique(vcat(envelope_indices(ch("F_tip"), 1400), envelope_indices(ch("F_root"), 1400),
            envelope_indices(ch("curvature_ratio_max"), 600))))
        channels = Dict{String,String}("t" => b64f32(t[idx]))
        for k in ("F_tip", "F_root", "tension_max", "curvature_ratio_max", "KE", "W_tip", "anchor_x", "progress")
            channels[k] = b64f32(ch(k)[idx])
        end

        tf = f["fields/t"]
        fi = collect(1:frame_stride:length(tf))
        X = f["fields/x"][fi, :]
        Z = f["fields/z"][fi, :]
        ni = collect(1:node_stride:size(X, 2))
        # positions in mm, laid out frame -> node -> (x, z) for JavaScript
        xy = zeros(Int16, 2, length(ni), length(fi))
        xy[1, :, :] = permutedims(round.(Int16, X[:, ni] .* 1000))
        xy[2, :, :] = permutedims(round.(Int16, Z[:, ni] .* 1000))
        # tension on each coarse segment: the larger-magnitude fine segment it covers
        T = f["fields/tension"][fi, :]
        Tc = zeros(Float32, length(ni) - 1, length(fi))
        for (k, _) in enumerate(fi), j in 1:length(ni)-1
            seg = view(T, k, ni[j]:ni[j+1]-1)
            Tc[j, k] = seg[argmax(abs.(seg))]
        end
        fold = permutedims(f["fields/fold_angle"][fi, :])
        hinge_nodes = f["layout/hinge_nodes"]

        Dict(
            "name" => name,
            "summary" => JSON.parse(f["summary_json"]),
            "channels" => channels,
            "frames" => Dict(
                "t" => b64f32(tf[fi]), "n" => length(fi), "nodes" => length(ni),
                "xy_mm" => b64i16(xy), "tension" => b64f32(Tc), "fold" => b64f32(fold),
                "n_hinges" => size(fold, 1),
                "hinge_idx" => [findfirst(==(h), ni) - 1 for h in hinge_nodes],
            ),
        )
    end
end

"""Read results/sweep/sweep_summary.csv into row dictionaries (numbers parsed)."""
function read_summary_csv(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    rows = Dict{String,Any}[]
    for line in lines[2:end]
        cells = String[]
        buf = IOBuffer()
        quoted = false
        for c in line
            if c == '"'
                quoted = !quoted
            elseif c == ',' && !quoted
                push!(cells, String(take!(buf)))
            else
                write(buf, c)
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

function main(cases)
    data = Dict{String,Any}("cases" => [export_case(c) for c in cases])
    sweep_csv = joinpath(ROOT, "results", "sweep", "sweep_summary.csv")
    isfile(sweep_csv) ? (data["sweep"] = read_summary_csv(sweep_csv)) : @warn "no sweep summary; the sensitivity section will be empty"
    template = read(joinpath(@__DIR__, "explorer_template.html"), String)
    occursin("<!--EXPLORER_DATA-->", template) || error("template is missing the <!--EXPLORER_DATA--> marker")
    # JSON cannot contain "</script" except inside strings; escape defensively.
    payload = replace(JSON.json(data), "</" => "<\\/")
    html = replace(template, "<!--EXPLORER_DATA-->" => "<script>window.FLAREWING_DATA = " * payload * ";</script>")
    out = joinpath(ROOT, "docs", "index.html")
    mkpath(dirname(out))
    write(out, html)
    @printf("wrote %s (%.1f MB, %d cases)\n", relpath(out, ROOT), filesize(out) / 1e6, length(cases))
end

main(isempty(ARGS) ? DEFAULT_CASES : ARGS)
