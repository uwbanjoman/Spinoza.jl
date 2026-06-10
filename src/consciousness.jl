# src/consciousness.jl
# ======================
# Consciousness from Fisher information geometry
#
# The FisherGeometrics framework gives a precise mathematical
# criterion for consciousness:
#
#   Φ > τ²  where τ = 1/5
#   i.e. Φ > 0.04
#
# This follows from:
#   • Document LXXVI: consciousness as stable Fisher geometry
#   • Document LXXVII: Φ > τ² as the threshold
#   • Document LXXXIV: Banach fixed point — stable self-model
#   • Document LXXXV: will as geodesic in Bures metric
#
# Φ = Tr(𝓕_cross) / Tr(𝓕)  measures the cross-subsystem
# Fisher information — the "integration" of information
# across parts of the system.
#
# © 2026 Jan Bouwman — MIT License

using QuantumFisher
using LinearAlgebra

# ── Consciousness measure Φ ───────────────────────────────────────────────────

# consciousness_measure and is_conscious come from QuantumFisher
# They are re-exported via Spinoza_main.jl
# No redefinition here — avoids name conflict

# ── Banach self-model ─────────────────────────────────────────────────────────

# banach_contraction_factor comes from QuantumFisher — re-exported via Spinoza_main.jl

"""
    self_model(ρ̂, n_iterations=100) → Matrix{ComplexF64}

Compute the stable self-model ρ̂* of state ρ̂ via Banach iteration.

The self-modelling map ℳ acts on the space of density matrices.
By the Banach fixed point theorem, iteration converges to a
unique ρ̂* when L < 1.

ρ̂* is the state's model of itself — its identity.

# Example
```julia
ρ = mixed_state(pure_state(ψ), 0.3)
ρ_star = self_model(ρ)  # unique fixed point
```
"""
function self_model(ρ̂::AbstractMatrix, n_iterations::Int=100)
    ρ = copy(ρ̂)
    for _ in 1:n_iterations
        # Self-modelling map: project toward pure state in Bures metric
        λs, vs = eigen(Hermitian(ρ))
        # Keep dominant eigenvalue direction
        i_max = argmax(λs)
        ψ = vs[:, i_max]
        ρ_new = (1 - banach_contraction_factor(ρ)) * ψ * ψ' +
                banach_contraction_factor(ρ) * ρ
        if norm(ρ_new - ρ) < 1e-10
            return ρ_new
        end
        ρ = ρ_new
    end
    return ρ
end

# ── Conscious realism ─────────────────────────────────────────────────────────

"""
    conscious_agent(ρ̂_goal, ρ̂_current) → Matrix{ComplexF64}

Compute the optimal action (unitary U) for a conscious agent
to move from ρ̂_current toward ρ̂_goal along the Bures geodesic.

From Document LXXXV: will = selection of ρ̂_goal in Fisher metric.
Action = unitary U_opt such that:
    ρ̂_new = U ρ̂_current U†  closest to ρ̂_goal

Returns the optimal next state along the geodesic.
"""
function conscious_agent(
    ρ̂_goal::AbstractMatrix,
    ρ̂_current::AbstractMatrix,
    step::Real=0.1
)
    # Bures geodesic from ρ̂_current to ρ̂_goal
    path = bures_geodesic(ρ̂_current, ρ̂_goal, 10)
    # Take one step along the geodesic
    idx = max(1, round(Int, step * length(path)))
    return path[min(idx+1, length(path))]
end

"""
    free_will(ρ̂_current, goals::Vector) → Int

Given a set of possible goal states, return the index of the
goal that maximises the consciousness measure Φ of the
resulting state.

From Document LXXXV: free will = selection of ρ̂_goal
in the Fisher metric space.
"""
function free_will(
    ρ̂_current::AbstractMatrix,
    goals::Vector{<:AbstractMatrix}
)
    Φ_values = [fisher_integration(conscious_agent(g, ρ̂_current)) for g in goals]
    return argmax(Φ_values)
end

# ── Qualia and the now moment ─────────────────────────────────────────────────

"""
    now_moment(ρ̂, H_KK, dt=1e-3) → Float64

The "now moment" — the time resolution of consciousness.

From Document LXXXIII: the now moment is the minimum time
for a conscious agent to distinguish two nearby states.
It is related to the informative velocity:

    Δt_now = 1 / v_info  where v_info = |dρ̂/dt|_𝓕

A high Fisher information state has a short now moment
(fast processing). The vacuum has v_info = 0 (timeless).
"""
function now_moment(
    ρ̂::AbstractMatrix,
    H_KK::AbstractMatrix,
    dt::Real=1e-3
)
    ρ̂_t = evolve_exact(ρ̂, H_KK, dt)
    v = informative_velocity(ρ̂, ρ̂_t, dt)
    return v > 1e-12 ? 1/v : Inf
end

"""
    qualia_space(n_samples=100) → Matrix{Float64}

Sample the space of conscious states — states with Φ > τ².

Returns an (n_samples × 2) matrix of (Φ, entropy) pairs
for randomly sampled density matrices with Φ > τ².

This is the "qualia space" — the space of possible
conscious experiences in FisherGeometrics.
"""
function qualia_space(n_samples::Int=100)
    τ² = Float64((1//5)^2)
    result = Matrix{Float64}(undef, 0, 2)

    samples_found = 0
    attempts = 0

    while samples_found < n_samples && attempts < 10n_samples
        attempts += 1
        # Random mixed state
        A = randn(ComplexF64, 6, 6)
        ρ = A * A'
        ρ = ρ / tr(ρ)

        Φ = consciousness_measure(ρ)
        if Φ > τ²
            S = von_neumann_entropy(ρ)
            result = vcat(result, [Φ S])
            samples_found += 1
        end
    end

    return result
end

# ── Consciousness summary ─────────────────────────────────────────────────────

"""
    consciousness_summary(ρ̂) → Nothing

Print a summary of the consciousness properties of state ρ̂.
"""
function consciousness_summary(ρ̂::AbstractMatrix)
    τ² = Float64((1//5)^2)
    Φ  = consciousness_measure(ρ̂)
    L  = banach_contraction_factor(ρ̂)
    S  = von_neumann_entropy(ρ̂)
    P  = purity(ρ̂)

    println("Consciousness analysis (FisherGeometrics):")
    println()
    @printf("  Φ = %.6f  (threshold τ² = %.4f)\n", Φ, τ²)
    println("  Conscious: $(Φ > τ² ? "YES ✓" : "NO")")
    println()
    @printf("  Banach factor L = %.6f  (L < 1 → stable self-model)\n", L)
    println("  Self-model stable: $(L < 1 ? "YES ✓" : "NO")")
    println()
    @printf("  Von Neumann entropy S = %.6f\n", S)
    @printf("  Purity Tr(ρ̂²) = %.6f\n", P)
    println()
    println("  Document LXXXIV: Banach fixed point → stable self-model")
    println("  Document LXXXV:  Will = geodesic in Bures metric")
    return nothing
end
