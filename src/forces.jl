# Force kernels. Each kernel is a pure function of its arguments: it ADDS
# nodal forces into caller-owned buffers and returns stored energy (or the
# force it applied). No kernel touches global or hidden state, so they can be
# reused for a 3D version (axial, attachment, gravity and settle damping are
# already dimension-generic over SVector{D}) and by parameter fitting.
#
# Units: forces N per metre width, energies J per metre width.
# Buffer convention: `F` receives the total force; `Fd` additionally receives
# the dissipative part so the integrator can book dissipation exactly.

"""
    axial_forces!(F, Fd, tension, x, v, rest_length, k, c) -> elastic energy

Spring-dashpot between consecutive nodes. `tension[j]` receives the total
axial force in segment j (spring + dashpot, positive in tension, N/m).
"""
function axial_forces!(F::AbstractVector{SVector{D,T}}, Fd, tension, x, v, rest_length, k, c) where {D,T}
    E = zero(T)
    @inbounds for j in eachindex(rest_length)
        d = x[j+1] - x[j]
        len = norm(d)
        e = d / len
        stretch = len - rest_length[j]
        rate = dot(e, v[j+1] - v[j])
        fs = k[j] * stretch
        fd = c[j] * rate
        f = (fs + fd) * e
        F[j] += f
        F[j+1] -= f
        Fd[j] += fd * e
        Fd[j+1] -= fd * e
        tension[j] = fs + fd
        E += 0.5k[j] * stretch^2
    end
    E
end

"""
    bending_forces!(F, Fd_panel, Fd_hinge, theta, curvature, x, v, lay,
                    k_bend, c_bend, law, fold_sign, hinge_state) -> (E_panel, E_hinge, E_released)

Torsional spring-dampers on the turning angle at interior nodes, converted to
nodal forces through the exact angle gradient. The angle is evaluated as
`atan2(cross, dot)` and unwrapped against the previous value stored in
`theta` (updated in place), so folds near +-pi are well defined and cannot
jump by 2pi. Panel nodes: `M = -k theta - c dtheta/dt` (flat rest shape).
Hinge nodes: the hinge law in fold coordinates `phi = fold_sign * theta`;
laws with internal state update `hinge_state` in place and report the energy
their update released (`E_released`, booked as dissipation).
`curvature[i]` receives theta / l_bar at panel-interior nodes (1/m), else 0.
"""
function bending_forces!(F, Fd_panel, Fd_hinge, theta, curvature, x, v, lay::Layout,
    k_bend, c_bend, law::HingeLaw, fold_sign, hinge_state)
    E_panel = 0.0
    E_hinge = 0.0
    E_released = 0.0
    N = lay.N
    @inbounds for i in 2:N-1
        kind = lay.kind[i]
        a = x[i] - x[i-1]
        b = x[i+1] - x[i]
        raw = atan(a[1] * b[2] - a[2] * b[1], a[1] * b[1] + a[2] * b[2])
        delta = raw - theta[i]
        if abs(delta) > pi                 # unwrap across the +-pi branch cut
            delta -= 2pi * round(delta / 2pi)
        end
        th = theta[i] + delta
        theta[i] = th
        # d(theta)/dx: perp(a)/|a|^2 at i-1, perp(b)/|b|^2 at i+1, minus both at i
        ga = SVector(-a[2], a[1]) / dot(a, a)
        gb = SVector(-b[2], b[1]) / dot(b, b)
        gim1 = ga
        gip1 = gb
        gi = -ga - gb
        thdot = dot(gim1, v[i-1]) + dot(gi, v[i]) + dot(gip1, v[i+1])
        if kind == PANEL_NODE
            Mel = -k_bend[i] * th
            Md = -c_bend[i] * thdot
            E_panel += 0.5k_bend[i] * th^2
            lbar = 0.5(lay.rest_length[i-1] + lay.rest_length[i])
            curvature[i] = th / lbar
            Mtot = Mel + Md
            F[i-1] += Mtot * gim1
            F[i] += Mtot * gi
            F[i+1] += Mtot * gip1
            Fd_panel[i-1] += Md * gim1
            Fd_panel[i] += Md * gi
            Fd_panel[i+1] += Md * gip1
        else # HINGE_NODE
            s = fold_sign[i]
            phi = s * th
            phidot = s * thdot
            st, released = hinge_update_state(law, phi, hinge_state[i])
            hinge_state[i] = st
            E_released += released
            Mphi = hinge_moment(law, phi, phidot, st)
            Md = -s * hinge_damping(law) * phidot   # damping part, theta coordinates
            Mtot = s * Mphi
            E_hinge += hinge_energy(law, phi, st)
            curvature[i] = 0.0
            F[i-1] += Mtot * gim1
            F[i] += Mtot * gi
            F[i+1] += Mtot * gip1
            Fd_hinge[i-1] += Md * gim1
            Fd_hinge[i] += Md * gi
            Fd_hinge[i+1] += Md * gip1
        end
    end
    (E_panel, E_hinge, E_released)
end

"""
    attachment_force(x, v, anchor, anchor_velocity, k, c) -> (F_spring, F_damper)

Zero-rest-length spring-damper from a node to an anchor point. Returns the
forces ON THE NODE (N/m); the load on the structure is their negative.
"""
@inline function attachment_force(x::SVector, v::SVector, a::SVector, adot::SVector, k, c)
    (-k * (x - a), -c * (v - adot))
end

attachment_energy(x::SVector, a::SVector, k) = 0.5k * dot(x - a, x - a)

"""
    gravity_forces!(F, mass, g, x) -> potential energy (relative to the origin)
"""
function gravity_forces!(F, mass, g::SVector, x)
    PE = 0.0
    @inbounds for i in eachindex(mass)
        F[i] += mass[i] * g
        PE -= mass[i] * dot(g, x[i])
    end
    PE
end

"""
    drag_forces!(F, Fd, x, v, tributary, rho, Cd)

Quadratic air drag on the velocity component normal to the local blanket
surface, `F = -0.5 rho Cd A v_n |v_n| n`, with A the node's tributary area per
metre of width. The surface normal at a node is taken from its neighbours.
"""
function drag_forces!(F, Fd, x::AbstractVector{SVector{2,T}}, v, tributary, rho, Cd) where {T}
    N = length(x)
    @inbounds for i in 1:N
        t = x[min(i + 1, N)] - x[max(i - 1, 1)]
        n = SVector(-t[2], t[1]) / norm(t)
        vn = dot(v[i], n)
        f = (-0.5 * rho * Cd * tributary[i] * vn * abs(vn)) * n
        F[i] += f
        Fd[i] += f
    end
    nothing
end

"""
    mass_damping_forces!(F, Fd, mass, v, alpha)

Mass-proportional damping `-alpha m v`, used only during the settle phase.
"""
function mass_damping_forces!(F, Fd, mass, v, alpha)
    @inbounds for i in eachindex(mass)
        f = -alpha * mass[i] * v[i]
        F[i] += f
        Fd[i] += f
    end
    nothing
end
