# Node layout (connectivity, rest lengths, lumped masses) and initial shapes.
#
# Coordinates: planar x-z. x is the deployment direction and the stacking
# direction of the stowed accordion; z is along the stowed panels and normal
# to the deployed blanket (gravity along -z sags the deployed blanket).
#
# Node numbering runs root (1) to tip (N). Panel k (1-based) spans nodes
# (k-1)*nps+1 .. k*nps+1; hinge k joins panels k and k+1 at node k*nps+1.

const END_NODE = Int8(0)
const PANEL_NODE = Int8(1)
const HINGE_NODE = Int8(2)

"""
    Layout

Connectivity and lumped properties of the discretized blanket. Masses and
tributary lengths are per metre of width.
"""
struct Layout
    n_panels::Int
    panel_length::Float64
    nodes_per_panel::Int          # segments per panel
    areal_density::Float64
    N::Int                        # number of nodes
    rest_length::Vector{Float64}  # per segment, m
    mass::Vector{Float64}         # per node, kg per m width
    tributary::Vector{Float64}    # per node, m
    kind::Vector{Int8}            # END_NODE | PANEL_NODE | HINGE_NODE
    panel::Vector{Int}            # owning panel (hinges report the root-side panel)
    hinge_nodes::Vector{Int}      # node index of hinge k (k = 1..n_panels-1)
    arc::Vector{Float64}          # unstretched arc length from root, m
end

function build_layout(n_panels::Int, panel_length::Float64, nps::Int, areal_density::Float64)
    N = n_panels * nps + 1
    l = panel_length / nps
    rest = fill(l, N - 1)
    trib = zeros(N)
    for j in 1:N-1
        trib[j] += 0.5rest[j]
        trib[j+1] += 0.5rest[j]
    end
    kind = fill(PANEL_NODE, N)
    kind[1] = END_NODE
    kind[N] = END_NODE
    hinges = [k * nps + 1 for k in 1:n_panels-1]
    kind[hinges] .= HINGE_NODE
    panel = [min(n_panels, max(1, cld(i - 1, nps))) for i in 1:N]
    panel[1] = 1
    arc = [(i - 1) * l for i in 1:N]
    Layout(n_panels, panel_length, nps, areal_density, N, rest, areal_density .* trib,
        trib, kind, panel, hinges, arc)
end

"""Indices of nodes strictly inside a panel (where cell curvature is checked)."""
panel_interior_nodes(lay::Layout) = findall(==(PANEL_NODE), lay.kind)

"""
    stowed_positions(lay, stow_gap) -> Vector{SVector{2}}

Z-fold stowed stack. Fold (hinge) points alternate between z = 0 and z = -h
and step `stow_gap` along x, so each straight panel leans by asin(gap/L).
Adjacent layers are on average `stow_gap` apart and every fold keeps a small
finite opening angle of 2*asin(gap/L): there is no zero-radius crease. Root
is at the origin. Hinge k has fold sign (-1)^(k+1) (see `SimParams`).
"""
function stowed_positions(lay::Layout, stow_gap::Float64)
    L = lay.panel_length
    stow_gap < L || error("stow_gap must be smaller than panel_length")
    h = sqrt(L^2 - stow_gap^2)
    nps = lay.nodes_per_panel
    x = Vector{SVector{2,Float64}}(undef, lay.N)
    corner(k) = SVector(k * stow_gap, iseven(k) ? 0.0 : -h)
    for k in 1:lay.n_panels
        a, b = corner(k - 1), corner(k)
        for j in 0:nps
            x[(k-1)*nps+1+j] = a + (j / nps) * (b - a)
        end
    end
    x
end

"""
    flat_positions(lay; strain=0, origin=(0,0), angle=0) -> Vector{SVector{2}}

Straight blanket along direction `angle` (radians from +x), uniformly
stretched by `strain`. Used by the verification tests.
"""
function flat_positions(lay::Layout; strain = 0.0, origin = SVector(0.0, 0.0), angle = 0.0)
    e = SVector(cos(angle), sin(angle))
    [SVector{2,Float64}(origin + (1 + strain) * lay.arc[i] * e) for i in 1:lay.N]
end

"""Signed turning angle at every interior node, atan2(cross, dot)."""
function turning_angles(x::AbstractVector{<:SVector{2}})
    N = length(x)
    th = zeros(N)
    for i in 2:N-1
        a = x[i] - x[i-1]
        b = x[i+1] - x[i]
        th[i] = atan(a[1] * b[2] - a[2] * b[1], a[1] * b[1] + a[2] * b[2])
    end
    th
end

"""Current length of the polyline through all nodes, m."""
polyline_length(x) = sum(norm(x[i+1] - x[i]) for i in 1:length(x)-1)
