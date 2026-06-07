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
using DifferentialEquations
using Printf

# Re-export from QuantumFisher
export vacuum_state, pure_state, mixed_state,
       von_neumann_entropy, purity,
       bures_distance, bures_fidelity, bures_geodesic,
       kk_hamiltonian, evolve_exact, informative_velocity,
       fisher_tensor, fisher_scalar,
       consciousness_measure, is_conscious, banach_contraction_factor

include("spacetime.jl")
include("smoothness.jl")
include("gravity.jl")
include("matter.jl")
include("consciousness.jl")
include("higgs.jl")
include("simulation.jl")

export
    # spacetime.jl
    spacetime_generators,
    covariance_metric,
    minkowski_metric,
    lorentz_interval,
    is_timelike,
    is_lightlike,
    is_spacelike,
 
    # smoothness.jl
    is_in_D6,
    D6_diameter,
    D6_bures_diameter,
    bgk_vn_rhs,
    lipschitz_constant,
    velocity_bound,
    gradient_bound,
    ns_smooth,
 
    # gravity.jl
    fisher_metric,
    metric_field,
    christoffel,
    ricci_tensor,
    ricci_scalar,
    einstein_tensor,
    stress_energy,
    einstein_equation_residual,
    gravitational_wave_metric,
 
    # matter.jl
    kk_spectrum_full,
    gauge_bosons,
    fermion_content,
    hypercharges,
    weinberg_angle,
    weinberg_angle_full,
    generation_count,
    particle_mass,
    cp_phase,
    cp_phase_degrees,
    sm_summary,
 
    # consciousness.jl
    self_model,
    conscious_agent,
    free_will,
    now_moment,
    qualia_space,
    consciousness_summary,

    # higgs.jl
    M_10D_planck,
    lyapunov_higgs,
 
    # simulation.jl
    simulate_1d,
    simulate_berry,
    simulate_3d,
    simulate_gravitational_waves,
    simulate_black_hole,
    simulate_cosmology

end # module Spinoza
