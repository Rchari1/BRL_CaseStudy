# Visual system for every figure: one palette, one type scale, recessive
# chrome. Colors come from a colorblind-validated reference palette (checked
# with a CVD/contrast validator): categorical slots are assigned in fixed
# order and follow the entity (tip is always slot 1, root always slot 2).

module Style

using CairoMakie
using CairoMakie: Makie

hex(s) = parse(Makie.Colors.Colorant, s)

# Surfaces and ink
const SURFACE = hex("#fcfcfb")
const PLANE = hex("#f9f9f7")
const INK = hex("#0b0b0b")
const INK2 = hex("#52514e")
const MUTED = hex("#898781")
const GRID = hex("#e1e0d9")
const BASELINE = hex("#c3c2b7")
const WASH = hex("#f0efec")

# Categorical slots (fixed order)
const BLUE = hex("#2a78d6")
const ORANGE = hex("#eb6834")
const AQUA = hex("#1baf7a")
const YELLOW = hex("#eda100")
const MAGENTA = hex("#e87ba4")
const GREEN = hex("#008300")
const VIOLET = hex("#4a3aa7")
const RED = hex("#e34948")
const SERIES = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]

# Entity colors used everywhere
const TIP = BLUE
const ROOT = ORANGE

# Status (reserved: only for pass/fail meaning, always with a label)
const CRITICAL = hex("#d03b3b")
const WARNING = hex("#fab219")
const GOOD = hex("#0ca30c")

# Sequential blue ramp (light -> dark)
const BLUE_RAMP = hex.(["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"])
const SEQUENTIAL = cgrad(vcat([SURFACE], BLUE_RAMP))

"""Ordinal ramp with n visible steps (validated light end, no near-surface step)."""
function ordinal(n::Int)
    steps = BLUE_RAMP[4:end]          # from #86b6ef (>= 2:1 on the surface) to #0d366b
    n == 1 && return [BLUE]
    [steps[round(Int, 1 + (length(steps) - 1) * (k - 1) / (n - 1))] for k in 1:n]
end

# Diverging: compression (blue) <- neutral gray -> tension (red)
const DIVERGING = cgrad(hex.(["#0d366b", "#1c5cab", "#3987e5", "#86b6ef", "#cde2fb",
    "#f0efec",
    "#fbd5d0", "#f4a39a", "#e34948", "#b52a2c", "#7a1416"]))

# Line variant: every step keeps contrast against the surface (slack = mid gray)
const DIVERGING_LINE = cgrad(hex.(["#0d366b", "#256abf", "#6da7ec", "#9a988f",
    "#ee8f80", "#e34948", "#b52a2c", "#7a1416"]))

const FONT_DIR = normpath(joinpath(@__DIR__, "..", "assets", "fonts"))

"""Bundled Inter (SIL OFL): identical rendering on every machine."""
function fonts()
    f(name) = joinpath(FONT_DIR, name)
    (regular = f("Inter-Regular.ttf"), medium = f("Inter-Medium.ttf"), bold = f("Inter-SemiBold.ttf"))
end

function theme()
    f = fonts()
    Theme(
        fonts = (regular = f.regular, bold = f.bold, medium = f.medium),
        fontsize = 13,
        backgroundcolor = SURFACE,
        textcolor = INK,
        linewidth = 2,
        Axis = (
            backgroundcolor = SURFACE,
            xgridcolor = GRID, ygridcolor = GRID, xgridwidth = 0.75, ygridwidth = 0.75,
            xgridvisible = false, ygridvisible = true,
            topspinevisible = false, rightspinevisible = false, leftspinevisible = false,
            bottomspinecolor = BASELINE, bottomspinewidth = 1,
            xtickcolor = BASELINE, ytickcolor = BASELINE, xticksize = 4, yticksize = 0,
            xticklabelcolor = INK2, yticklabelcolor = INK2,
            xticklabelsize = 11.5, yticklabelsize = 11.5,
            xlabelcolor = INK2, ylabelcolor = INK2, xlabelsize = 12, ylabelsize = 12,
            xlabelpadding = 6, ylabelpadding = 8,
            titlefont = :medium, titlesize = 13.5, titlecolor = INK, titlealign = :left, titlegap = 8,
            subtitlecolor = INK2, subtitlesize = 11.5, subtitlegap = 2,
            xminorticksvisible = false, yminorticksvisible = false,
        ),
        Colorbar = (
            spinewidth = 0, ticklabelcolor = INK2, ticklabelsize = 11, labelcolor = INK2, labelsize = 12,
            tickcolor = BASELINE, ticksize = 3, size = 10,
        ),
        Legend = (
            framevisible = false, labelcolor = INK2, labelsize = 12, patchsize = (16, 10),
            rowgap = 2, colgap = 14, padding = (0, 0, 0, 0), backgroundcolor = :transparent,
        ),
        Label = (color = INK,),
        Lines = (linewidth = 2, joinstyle = :round, linecap = :round),
    )
end

end # module Style
