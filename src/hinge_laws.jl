# Fold-hinge constitutive laws: THE HOOK FOR TEST DATA.
#
# A hinge law gives the moment acting on a fold as a function of its fold angle
# phi (rad; pi = fully stowed, 0 = flat, same sign convention at every hinge)
# and fold rate. Moments are per metre of width, positive when they act to
# INCREASE phi. To use measured data, either point `hinge.table_file` at a
# moment-angle CSV (TabulatedHinge) or add a new subtype of `HingeLaw` and
# define `hinge_moment` and `hinge_energy` for it. Hysteretic or
# viscoelastic laws with internal state belong here too (stretch goal).

abstract type HingeLaw end

"""
    LinearHinge(k, phi_rest, c)

`M = -k (phi - phi_rest) - c dphi/dt`. `phi_rest` between 0 (flat) and pi
(stowed) represents fold memory.
"""
struct LinearHinge <: HingeLaw
    k::Float64
    phi_rest::Float64
    c::Float64
end

hinge_moment(h::LinearHinge, phi::Float64, phidot::Float64) = -h.k * (phi - h.phi_rest) - h.c * phidot
hinge_energy(h::LinearHinge, phi::Float64) = 0.5h.k * (phi - h.phi_rest)^2
hinge_stiffness_bound(h::LinearHinge) = h.k
hinge_damping(h::LinearHinge) = h.c

"""
    TabulatedHinge(phi, M, c)
    TabulatedHinge(csv_path, c)

Piecewise-linear elastic moment-angle curve from test data plus linear
damping. `M[i]` is the elastic moment at fold angle `phi[i]` (sorted
ascending); it is extrapolated linearly outside the table. The CSV has two
columns, `fold_angle_rad,moment_Nm_per_m`, with an optional header row.
Stored energy is the exact integral of the curve (zero at the first point).
"""
struct TabulatedHinge <: HingeLaw
    phi::Vector{Float64}
    M::Vector{Float64}
    E::Vector{Float64}   # -integral of M from phi[1] to phi[i]
    c::Float64
    kmax::Float64
end

function TabulatedHinge(phi::Vector{Float64}, M::Vector{Float64}, c::Float64)
    length(phi) == length(M) >= 2 || error("hinge table needs at least two rows")
    issorted(phi) && allunique(phi) || error("hinge table angles must be strictly increasing")
    E = zeros(length(phi))
    for i in 2:length(phi)
        E[i] = E[i-1] - 0.5(M[i] + M[i-1]) * (phi[i] - phi[i-1])
    end
    kmax = maximum(abs.(diff(M) ./ diff(phi)))
    TabulatedHinge(phi, M, E, c, kmax)
end

function TabulatedHinge(path::AbstractString, c::Float64)
    rows = readdlm(path, ',')
    first_numeric = findfirst(r -> rows[r, 1] isa Real, 1:size(rows, 1))
    data = Float64.(rows[first_numeric:end, 1:2])
    TabulatedHinge(data[:, 1], data[:, 2], c)
end

function _segment(h::TabulatedHinge, phi::Float64)
    n = length(h.phi)
    clamp(searchsortedlast(h.phi, phi), 1, n - 1)
end

function hinge_moment(h::TabulatedHinge, phi::Float64, phidot::Float64)
    i = _segment(h, phi)
    w = (phi - h.phi[i]) / (h.phi[i+1] - h.phi[i])
    (1 - w) * h.M[i] + w * h.M[i+1] - h.c * phidot
end

function hinge_energy(h::TabulatedHinge, phi::Float64)
    i = _segment(h, phi)
    slope = (h.M[i+1] - h.M[i]) / (h.phi[i+1] - h.phi[i])
    d = phi - h.phi[i]
    h.E[i] - (h.M[i] * d + 0.5slope * d^2)
end

hinge_stiffness_bound(h::TabulatedHinge) = h.kmax
hinge_damping(h::TabulatedHinge) = h.c
