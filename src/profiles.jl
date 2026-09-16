# Prescribed tip-attachment motion. A profile maps time since deployment start
# to normalized progress s in [0, 1] and its rate ds/dt. The anchor moves on a
# straight line from its stowed position to its deployed end point.

abstract type DeployProfile end

"""Quintic S-curve: zero velocity and acceleration at both ends."""
struct SmoothProfile <: DeployProfile
    T::Float64
end

"""Constant acceleration for `f*T`, constant speed, constant deceleration for `f*T`."""
struct TrapezoidProfile <: DeployProfile
    T::Float64
    f::Float64
end

"""
Smooth velocity ramp over `f*T`, then constant speed and an abrupt stop at the
end position (worst-case snap-taut). `f = 0` also makes the start abrupt.
"""
struct ConstantThenStopProfile <: DeployProfile
    T::Float64
    f::Float64
end

profile_name(::SmoothProfile) = "smooth"
profile_name(::TrapezoidProfile) = "trapezoid"
profile_name(::ConstantThenStopProfile) = "constant_then_stop"

function make_profile(name::AbstractString, T::Float64, accel_fraction::Float64, ramp_fraction::Float64)
    name = lowercase(name)
    if name == "smooth"
        SmoothProfile(T)
    elseif name == "trapezoid"
        0 < accel_fraction <= 0.5 || error("deployment.accel_fraction must be in (0, 0.5]")
        TrapezoidProfile(T, accel_fraction)
    elseif name == "constant_then_stop"
        0 <= ramp_fraction < 1 || error("deployment.start_ramp_fraction must be in [0, 1)")
        ConstantThenStopProfile(T, ramp_fraction)
    else
        error("deployment.profile must be smooth, trapezoid or constant_then_stop (got `$name`)")
    end
end

"""
    progress(profile, t) -> (s, sdot)

Normalized progress and its time derivative at time `t` (s) after deployment
start. Holds s = 0 for t < 0 and s = 1 after the profile ends.
"""
function progress(p::SmoothProfile, t::Float64)
    t <= 0 && return (0.0, 0.0)
    t >= p.T && return (1.0, 0.0)
    u = t / p.T
    s = u^3 * (10 - 15u + 6u^2)
    sdot = 30u^2 * (1 - u)^2 / p.T
    (s, sdot)
end

function progress(p::TrapezoidProfile, t::Float64)
    t <= 0 && return (0.0, 0.0)
    t >= p.T && return (1.0, 0.0)
    u, f = t / p.T, p.f
    vmax = 1 / (1 - f)                       # in units of 1/T
    if u < f
        (0.5vmax * u^2 / f, vmax * u / f / p.T)
    elseif u <= 1 - f
        (vmax * (u - 0.5f), vmax / p.T)
    else
        r = 1 - u
        (1 - 0.5vmax * r^2 / f, vmax * r / f / p.T)
    end
end

function progress(p::ConstantThenStopProfile, t::Float64)
    t <= 0 && return (0.0, 0.0)
    t >= p.T && return (1.0, 0.0)
    u, f = t / p.T, p.f
    V = 1 / (1 - 0.5f)
    if f > 0 && u < f
        r = u / f
        (V * f * (r^3 - 0.5r^4), V * (3r^2 - 2r^3) / p.T)
    else
        (V * (0.5f + (u - f)), V / p.T)
    end
end

"""
    AnchorMotion

Straight-line anchor path from `a0` to `a1` following `profile`.
`anchor_state(m, t)` returns position and velocity.
"""
struct AnchorMotion{P<:DeployProfile}
    a0::SVector{2,Float64}
    a1::SVector{2,Float64}
    profile::P
end

function anchor_state(m::AnchorMotion, t::Float64)
    s, sdot = progress(m.profile, t)
    d = m.a1 - m.a0
    (m.a0 + s * d, sdot * d)
end

"""Anchor held at a fixed point."""
struct FixedAnchor
    a::SVector{2,Float64}
end
anchor_state(m::FixedAnchor, ::Float64) = (m.a, SVector(0.0, 0.0))

"""Anchor driven by an arbitrary function `t -> (position, velocity)` (tests, future coupling)."""
struct FunctionAnchor{F}
    f::F
end
anchor_state(m::FunctionAnchor, t::Float64) = m.f(t)

"""
    deployed_tip_position(p)

End point of the tip anchor: flat blanket length x (1 + preload_strain) along
+x from the root anchor.
"""
deployed_tip_position(p) = SVector(flat_length(p) * (1 + p.preload_strain), 0.0)
