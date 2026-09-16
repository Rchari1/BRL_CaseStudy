# Fold-hinge constitutive laws: THE HOOK FOR TEST DATA.
#
# A hinge law gives the moment acting on a fold as a function of its fold angle
# phi (rad; pi = fully stowed, 0 = flat, same sign convention at every hinge)
# and fold rate. Moments are per metre of width, positive when they act to
# INCREASE phi. To use measured data, either point `hinge.table_file` at a
# moment-angle CSV (TabulatedHinge) or add a new subtype of `HingeLaw` and
# define `hinge_moment` and `hinge_energy` for it. Laws with internal state
# (hysteresis, creep) also define `hinge_state_init` and `hinge_update_state`;
# see `HystereticHinge`.

abstract type HingeLaw end

# Stateless laws ignore the per-hinge state argument.
hinge_moment(h::HingeLaw, phi::Float64, phidot::Float64, state::Float64) = hinge_moment(h, phi, phidot)
hinge_energy(h::HingeLaw, phi::Float64, state::Float64) = hinge_energy(h, phi)
"""Initial per-hinge state for a fold starting at angle `phi0` (stateless laws: unused)."""
hinge_state_init(::HingeLaw, phi0::Float64) = 0.0
"""
    hinge_update_state(law, phi, state) -> (new_state, dissipated_energy)

Advance internal state at fold angle `phi`. Returns the energy (J per metre
width) released irreversibly by the update, which is booked as dissipation.
"""
hinge_update_state(::HingeLaw, phi::Float64, state::Float64) = (state, 0.0)

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

"""
    HystereticHinge(k_rest, phi_rest, k_elastic, yield_moment, c)

Jenkins-type fold with hysteresis (stretch goal: stand-in for viscoelastic fold
memory). A spring `k_rest` to the rest angle (fold memory) acts in parallel with
a spring `k_elastic` in series with a Coulomb slider that slips at
`yield_moment`. Small reversals are stiff (`k_rest + k_elastic`); once the
slider slips the fold follows `k_rest`, so loading and unloading take
different paths and the loop area is dissipated. The slider angle is the
per-hinge state and starts at the stowed fold angle (the fold is "set" by
stowage). Energy released when the slider slips is booked as hinge dissipation.
"""
struct HystereticHinge <: HingeLaw
    k_rest::Float64
    phi_rest::Float64
    k_elastic::Float64
    yield_moment::Float64
    c::Float64
end

hinge_moment(h::HystereticHinge, phi::Float64, phidot::Float64, phi_s::Float64) =
    -h.k_rest * (phi - h.phi_rest) - h.k_elastic * (phi - phi_s) - h.c * phidot
hinge_energy(h::HystereticHinge, phi::Float64, phi_s::Float64) =
    0.5h.k_rest * (phi - h.phi_rest)^2 + 0.5h.k_elastic * (phi - phi_s)^2
hinge_state_init(::HystereticHinge, phi0::Float64) = phi0
function hinge_update_state(h::HystereticHinge, phi::Float64, phi_s::Float64)
    h.k_elastic > 0 || return (phi_s, 0.0)
    slip = h.yield_moment / h.k_elastic
    d = phi - phi_s
    new = d > slip ? phi - slip : d < -slip ? phi + slip : phi_s
    released = 0.5h.k_elastic * (d^2 - (phi - new)^2)
    (new, released)
end
hinge_moment(h::HystereticHinge, phi::Float64, phidot::Float64) = error("HystereticHinge needs its state; call hinge_moment(law, phi, phidot, state)")
hinge_energy(h::HystereticHinge, phi::Float64) = error("HystereticHinge needs its state")
hinge_stiffness_bound(h::HystereticHinge) = h.k_rest + h.k_elastic
hinge_damping(h::HystereticHinge) = h.c
