#!/usr/bin/env julia
# Command-line entry point.
#
#   julia --project=. -t auto bin/blanket_sim.jl run   configs/baseline.yaml
#   julia --project=. -t auto bin/blanket_sim.jl sweep configs/sweep.yaml
#
# Run from the repository root. All results are placeholder-parameter estimates.

using BlanketSim
exit(BlanketSim.main(ARGS))
