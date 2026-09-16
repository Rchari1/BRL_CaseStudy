# Verification tests against analytic cases (build brief, Section 6).
# These check the numerics, not the placeholder physics.

using Test
using LinearAlgebra
using StaticArrays
using BlanketSim

"""Config dictionary from the defaults plus dotted-key overrides."""
function cfg_with(pairs::Pair...)
    cfg = default_config()
    cfg["run_name"] = "test"
    for (k, v) in pairs
        set_param!(cfg, k, v)
    end
    cfg
end

"""Times at which `y` crosses `level` upward, linearly interpolated."""
function upward_crossings(t, y, level = 0.0)
    out = Float64[]
    for k in 1:length(y)-1
        if y[k] < level <= y[k+1]
            push!(out, t[k] + (level - y[k]) / (y[k+1] - y[k]) * (t[k+1] - t[k]))
        end
    end
    out
end

"""First time `y` reaches `level`, linearly interpolated."""
function first_crossing(t, y, level)
    k = findfirst(>=(level), y)
    k === nothing && return NaN
    k == 1 && return t[1]
    t[k-1] + (level - y[k-1]) / (y[k] - y[k-1]) * (t[k] - t[k-1])
end

rel(a, b) = abs(a - b) / abs(b)

const QUIET = (verbose = false,)

@testset "BlanketSim verification" begin

    @testset "1. Compound pendulum period" begin
        L, g, amp = 0.5, 9.80665, deg2rad(2.0)
        cfg = cfg_with("blanket.n_panels" => 1, "blanket.panel_length" => L,
            "blanket.nodes_per_panel" => 20, "blanket.panel_EI" => 10.0,
            "attachments.k_att_root" => 1e6, "attachments.c_att_root" => 0.0,
            "attachments.k_att_tip" => 0.0, "attachments.c_att_tip" => 0.0,
            "environment.gravity" => true, "environment.gravity_direction" => "-z",
            "deployment.T_settle" => 0.0, "output.field_dt" => 1e-3)
        p = SimParams(cfg)
        x0 = flat_positions(p.layout; angle = -pi / 2 + amp)
        r = simulate(p; x0, tip_motion = FixedAnchor(x0[end]), t_start = 0.0, t_end = 6.0, QUIET...)
        tx = r.ch["tip_x"]
        crossings = upward_crossings(r.t_ch, tx .- 0.0)
        period = (crossings[end] - crossings[1]) / (length(crossings) - 1)
        expected = 2pi * sqrt(2L / (3g))       # uniform rod pivoting about one end
        @info "pendulum" period expected error = rel(period, expected)
        @test length(crossings) >= 4
        @test rel(period, expected) < 0.01
    end

    @testset "2. Hinge oscillator frequency" begin
        # Two identical rigid panels joined by one hinge, floating free, gravity
        # off. With zero total momentum and angular momentum the fold angle
        # obeys I phi'' = -k phi with I = rho L^3 / 24 (each panel contributes
        # its centroidal inertia rho L^3/12 and the pair shares one rotation).
        L, rho, k, phi0 = 0.5, 1.0, 2.6e-2, 0.05
        cfg = cfg_with("blanket.n_panels" => 2, "blanket.panel_length" => L,
            "blanket.nodes_per_panel" => 20, "blanket.panel_EI" => 10.0,
            "blanket.areal_density" => rho,
            "hinge.k_hinge" => k, "hinge.theta_rest" => 0.0, "hinge.damping_ratio" => 0.0,
            "attachments.k_att_root" => 0.0, "attachments.c_att_root" => 0.0,
            "attachments.k_att_tip" => 0.0, "attachments.c_att_tip" => 0.0,
            "deployment.T_settle" => 0.0, "output.field_dt" => 1e-3)
        p = SimParams(cfg)
        lay = p.layout
        nps = lay.nodes_per_panel
        x0 = flat_positions(lay)
        hinge = x0[nps+1]
        for i in nps+2:lay.N
            d = x0[i] - hinge
            x0[i] = hinge + SVector(cos(phi0) * d[1], sin(phi0) * d[1])
        end
        r = simulate(p; x0, t_start = 0.0, t_end = 8.0, QUIET...)
        phi = r.fold[:, 1]
        crossings = upward_crossings(r.t_f, phi)
        period = (crossings[end] - crossings[1]) / (length(crossings) - 1)
        omega = 2pi / period
        expected = sqrt(k / (rho * L^3 / 24))
        @info "hinge oscillator" omega expected error = rel(omega, expected)
        @test length(crossings) >= 3
        @test rel(omega, expected) < 0.01
    end

    @testset "3. Energy conservation (no damping, no prescribed motion)" begin
        cfg = cfg_with("blanket.n_panels" => 4, "blanket.axial_damping_ratio" => 0.0,
            "blanket.panel_bending_damping_ratio" => 0.0, "hinge.damping_ratio" => 0.0,
            "attachments.c_att_root" => 0.0, "attachments.c_att_tip" => 0.0,
            "deployment.T_settle" => 0.0, "deployment.settle_mass_damping" => 0.0)
        p = SimParams(cfg)
        x0 = stowed_positions(p.layout, p.stow_gap)
        r = simulate(p; x0, tip_motion = FixedAnchor(x0[end]), t_start = 0.0, t_end = 4.0, QUIET...)
        E = r.ch["KE"] .+ r.ch["E_axial"] .+ r.ch["E_panel"] .+ r.ch["E_hinge"] .+ r.ch["E_att"] .+ r.ch["PE_grav"]
        drift = maximum(abs.(E .- E[1])) / E[1]
        moved = maximum(r.ch["KE"]) / E[1]
        @info "energy conservation" drift E0 = E[1] max_KE_fraction = moved
        @test moved > 0.05          # the hinges really do drive motion
        @test drift < 0.005
    end

    @testset "4. Energy balance with prescribed motion" begin
        for profile in ("smooth", "constant_then_stop")
            cfg = cfg_with("blanket.n_panels" => 4, "deployment.profile" => profile,
                "deployment.T_deploy" => 3.0, "deployment.T_settle" => 0.5, "deployment.T_hold" => 1.0)
            p = SimParams(cfg)
            r = simulate(p; QUIET...)
            W = r.ch["W_tip"][end]
            D = sum(r.ch[k][end] for k in ("D_axial", "D_panel", "D_hinge", "D_att", "D_drag", "D_settle"))
            @info "energy balance" profile W D error = r.balance_error
            @test W > 0
            @test r.balance_error < 0.01
        end
    end

    @testset "5. Static tension = EA x strain" begin
        L_total = 1.0
        delta = 2e-3
        # Flat rest angle: this checks the axial law alone (fold memory would kink the blanket).
        cfg = cfg_with("blanket.n_panels" => 2, "blanket.panel_length" => L_total / 2,
            "hinge.theta_rest" => 0.0, "deployment.T_settle" => 2.0, "deployment.settle_mass_damping" => 20.0,
            "output.field_dt" => 0.05)
        p = SimParams(cfg)
        x0 = flat_positions(p.layout)
        tip = FixedAnchor(SVector(L_total + delta, 0.0))
        r = simulate(p; x0, tip_motion = tip, t_start = -2.0, t_end = 0.0, QUIET...)
        xs = r.pos[end, :, 1]
        strain = (xs[end] - xs[1]) / L_total - 1
        T = r.tension[end, :]
        expected_series = delta / (1 / p.k_root + L_total / p.EA + 1 / p.k_tip)
        @info "static tension" strain mean_T = sum(T) / length(T) EA_strain = p.EA * strain expected_series
        @test all(rel.(T, p.EA * strain) .< 1e-3)
        @test rel(sum(T) / length(T), expected_series) < 1e-3
    end

    @testset "6. Axial wave speed" begin
        rho, EA, v0 = 1.0, 1.25e5, 0.01
        cfg = cfg_with("blanket.n_panels" => 4, "blanket.panel_length" => 0.5,
            "blanket.nodes_per_panel" => 40, "blanket.axial_damping_ratio" => 0.0,
            "attachments.k_att_root" => 1e7, "attachments.c_att_root" => 0.0,
            "attachments.k_att_tip" => 1e7, "attachments.c_att_tip" => 0.0,
            "deployment.T_settle" => 0.0, "output.field_dt" => 2e-5, "output.channel_dt" => 1e-5)
        p = SimParams(cfg)
        lay = p.layout
        x0 = flat_positions(lay)
        a0 = x0[end]
        tip = FunctionAnchor(t -> t <= 0 ? (a0, SVector(0.0, 0.0)) : (a0 + SVector(v0 * t, 0.0), SVector(v0, 0.0)))
        c = sqrt(EA / rho)
        L = 2.0
        r = simulate(p; x0, tip_motion = tip, t_start = 0.0, t_end = 0.9L / c, QUIET...)
        T_incident = rho * c * v0
        mid = 0.5 * (lay.arc[1:end-1] + lay.arc[2:end])
        j_near = argmin(abs.(mid .- 1.5))     # 0.5 m from the tip
        j_far = argmin(abs.(mid .- 0.5))      # 1.5 m from the tip
        t_near = first_crossing(r.t_f, r.tension[:, j_near], 0.5T_incident)
        t_far = first_crossing(r.t_f, r.tension[:, j_far], 0.5T_incident)
        speed = (mid[j_near] - mid[j_far]) / (t_far - t_near)
        # Impedance check: pulling at v0 the time-averaged tip load is rho*c*v0
        # (impulse = momentum of the moving region). A pointwise plateau would be
        # biased by the ringing a discrete chain shows right behind a step.
        mean_tip_load = sum(r.ch["F_tip"]) / length(r.ch["F_tip"])
        @info "wave speed" speed expected = c error = rel(speed, c) mean_tip_load T_incident
        @test rel(speed, c) < 0.05
        @test rel(mean_tip_load, T_incident) < 0.05
    end

    @testset "7. Convergence of peak attachment force" begin
        # A short, fast deployment whose peak is set by the snap-taut event itself.
        # (Peaks in long, lightly damped runs also depend on the phase of residual
        # flapping; README "Numerical notes" covers that separately.)
        function peak(profile; nps = 20, dt_factor = 1.0)
            cfg = cfg_with("blanket.n_panels" => 2, "blanket.nodes_per_panel" => nps,
                "deployment.profile" => profile, "deployment.T_deploy" => 1.0,
                "deployment.T_settle" => 0.5, "deployment.T_hold" => 0.3, "output.field_dt" => 0.05)
            p0 = SimParams(cfg)
            set_param!(cfg, "numerics.dt", p0.dt * dt_factor)
            r = simulate(SimParams(cfg); record_fields = false, QUIET...)
            max(r.peaks.F_tip, r.peaks.F_root)
        end
        for profile in ("smooth", "constant_then_stop")
            base = peak(profile)
            finer_mesh = peak(profile; nps = 40)
            finer_dt = peak(profile; dt_factor = 0.5)
            @info "convergence" profile base mesh_change = rel(finer_mesh, base) dt_change = rel(finer_dt, base)
            @test rel(finer_mesh, base) < 0.02
            @test rel(finer_dt, base) < 0.02
        end
    end

    @testset "8. Geometry sanity" begin
        p = SimParams(default_config())
        lay = p.layout
        x0 = stowed_positions(lay, p.stow_gap)
        @test polyline_length(x0) ≈ lay.n_panels * lay.panel_length rtol = 1e-12
        th = turning_angles(x0)
        signs = [sign(th[i]) for i in lay.hinge_nodes]
        @test all(signs[k] == -signs[k+1] for k in 1:length(signs)-1)
        @test all(abs(th[i]) > 0.99pi for i in lay.hinge_nodes)
        @test all(abs(th[i]) < 1e-9 for i in 2:lay.N-1 if !(i in lay.hinge_nodes))
        @test all(sign(th[i]) == p.fold_sign[i] for i in lay.hinge_nodes)
        @test x0[end][1] ≈ lay.n_panels * p.stow_gap
        @test x0[end][2] ≈ 0.0 atol = 1e-12        # even panel count: tip on the root line

        # Fully deployed: anchor end point and the tip node after a short deployment.
        cfg = default_config()
        set_param!(cfg, "blanket.n_panels", 4)
        set_param!(cfg, "deployment.T_deploy", 4.0)
        set_param!(cfg, "deployment.T_settle", 1.0)
        set_param!(cfg, "deployment.T_hold", 2.0)
        p4 = SimParams(cfg)
        r = simulate(p4; verbose = false)
        L_end = flat_length(p4) * (1 + p4.preload_strain)
        @test r.ch["anchor_x"][end] ≈ L_end rtol = 1e-12
        tip_offset = r.ch["F_tip"][end] / p4.k_tip
        @test abs(r.ch["tip_x"][end] - L_end) < tip_offset + 1e-4
        # Deployed arc length = flat length stretched by the mean axial strain.
        arc_end = polyline_length([SVector(r.pos[end, i, 1], r.pos[end, i, 2]) for i in 1:p4.layout.N])
        mean_strain = sum(r.tension[end, :]) / size(r.tension, 2) / p4.EA
        @test arc_end ≈ flat_length(p4) * (1 + mean_strain) rtol = 2e-4
        @test arc_end >= L_end - tip_offset
    end

    @testset "Hysteretic hinge (stretch goal)" begin
        # Closed cycle of the fold angle: work done on the fold equals the energy
        # the law reports as released, and the loop encloses positive area.
        h = HystereticHinge(2.6e-2, 0.3pi, 0.2, 0.01, 0.0)
        A, n = 0.4, 4000
        path = vcat(range(0, A; length = n), range(A, -A; length = 2n), range(-A, A; length = 2n))
        s_ = 0.0
        work = 0.0
        released = 0.0
        E_start = NaN
        started = false
        for k in 2:length(path)
            φ0, φ1 = path[k-1], path[k]
            if !started && k == n + 1                 # steady loop starts at phi = A
                started = true
                E_start = hinge_energy(h, φ0, s_)
            end
            M0 = hinge_moment(h, φ0, 0.0, s_)
            s_, rel = hinge_update_state(h, φ1, s_)
            M1 = hinge_moment(h, φ1, 0.0, s_)
            if started
                work += -0.5(M0 + M1) * (φ1 - φ0)     # work done ON the fold
                released += rel
            end
        end
        loop_energy = hinge_energy(h, path[end], s_) - E_start
        @info "hysteretic hinge" work released loop_energy
        @test released > 0
        @test abs(loop_energy) < 1e-6
        @test isapprox(work, released; rtol = 2e-3)

        cfg = cfg_with("blanket.n_panels" => 4, "deployment.T_deploy" => 3.0, "deployment.T_settle" => 0.5,
            "deployment.T_hold" => 1.0, "hinge.law" => "hysteretic", "hinge.k_elastic" => 0.2,
            "hinge.yield_moment" => 0.01)
        r = simulate(SimParams(cfg); record_fields = false, QUIET...)
        cfg_lin = cfg_with("blanket.n_panels" => 4, "deployment.T_deploy" => 3.0, "deployment.T_settle" => 0.5,
            "deployment.T_hold" => 1.0)
        r_lin = simulate(SimParams(cfg_lin); record_fields = false, QUIET...)
        @info "hysteretic deployment" balance = r.balance_error D_hinge = r.ch["D_hinge"][end] D_hinge_linear = r_lin.ch["D_hinge"][end]
        @test r.balance_error < 0.01
        @test r.ch["D_hinge"][end] > r_lin.ch["D_hinge"][end]
    end

    @testset "Profiles and hinge-law hook" begin
        for prof in (SmoothProfile(10.0), TrapezoidProfile(10.0, 0.2), ConstantThenStopProfile(10.0, 0.1))
            @test progress(prof, -1.0) == (0.0, 0.0)
            @test progress(prof, 10.0)[1] ≈ 1.0
            @test progress(prof, 9.999999)[1] ≈ 1.0 atol = 1e-5
            h = 1e-6
            for t in (0.5, 3.0, 7.7)
                s1, _ = progress(prof, t + h)
                s0, _ = progress(prof, t - h)
                @test (s1 - s0) / 2h ≈ progress(prof, t)[2] rtol = 1e-4
            end
        end
        lin = LinearHinge(2.6e-2, 0.3pi, 0.0)
        phis = collect(range(0, pi, length = 5))
        tab = TabulatedHinge(phis, [hinge_moment(lin, φ, 0.0) for φ in phis], 0.0)
        for φ in (0.1, 1.3, 2.9)
            @test hinge_moment(tab, φ, 0.0) ≈ hinge_moment(lin, φ, 0.0)
            @test hinge_energy(tab, φ) - hinge_energy(tab, 0.3pi) ≈ hinge_energy(lin, φ) atol = 1e-12
        end
    end
end
