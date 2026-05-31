# ===============
# Spinoza.jl
#
# One density matrix. One postulate. Everything.
#
# A unified framework deriving the Standard Model,
# spacetime, gravity, and consciousness from the
# single postulate g_AB = 𝓕_AB/ρ₀.
#
# © 2026 Jan Bouwman — MIT License

module Spinoza

using QuantumFisher
using LinearAlgebra

# Re-export from QuantumFisher
using QuantumFisher: density_matrix, vacuum_state, fisher_tensor,
                     bures_distance, evolve_exact, kk_hamiltonian

include("spacetime.jl")
# include("gravity.jl")       # coming: Einstein tensor, curvature
# include("matter.jl")        # coming: KK modes, SM particles
# include("consciousness.jl") # coming: Hoffman, Banach, qualia
# include("simulation.jl")    # coming: lattice, gravitational waves

export
    # spacetime.jl
    spacetime_generators,
    covariance_metric,
    minkowski_metric,
    lorentz_interval,
    is_timelike,
    is_lightlike,
    is_spacelike

end # module Spinoza
