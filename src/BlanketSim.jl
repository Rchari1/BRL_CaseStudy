"""
    BlanketSim

First-fidelity planar (x-z) lumped-mass simulation of a z-folded solar blanket
pulled open by a prescribed tip motion. Predicts attachment loads, blanket
tension, fold angles and cell curvature.

ALL MATERIAL PARAMETERS ARE PLACEHOLDERS. Results show how loads scale; they
are not validated against test data.
"""
module BlanketSim

using LinearAlgebra
using Printf
using Statistics
using DelimitedFiles
using StaticArrays
using YAML
using JSON
using JLD2

include("hinge_laws.jl")
include("geometry.jl")
include("profiles.jl")
include("config.jl")
include("forces.jl")
include("integrator.jl")
include("outputs.jl")
include("theme.jl")
include("sweep.jl")
include("plots.jl")
include("cli.jl")

export default_config, load_config, set_param!, get_param, SimParams, print_startup
export Layout, build_layout, stowed_positions, flat_positions, turning_angles, polyline_length
export SmoothProfile, TrapezoidProfile, ConstantThenStopProfile, progress, AnchorMotion,
    FixedAnchor, FunctionAnchor, anchor_state, deployed_tip_position
export LinearHinge, TabulatedHinge, HystereticHinge, hinge_moment, hinge_energy, hinge_update_state
export simulate, SimResult, summary_metrics, summary_line, save_results, write_summary
export flat_length, blanket_mass
export plot_run, animate_run
export load_sweep, run_sweep, aggregate, sensitivity_table, plot_sweep, cli_run, cli_sweep

end
