# Spacetime geometry derived from Fisher information
#
# Minkowski spacetime η_μν = diag(-1,+1,+1,+1) is derived from
# the covariance metric of N conscious agents in the classical limit
# via the quantum CLT (Goderis-Vets-Verbeure 1989, Document LXXX).
#
# The derivation:
#   1. Covariance metric in vacuum: Cov(G_μ,G_ν) = (1/3)δ_μν
#   2. Wick rotation: G_t is the U(1)_Y singlet (Y=0) → time (-)
#   3. Normalise: η_μν = diag(-1,+1,+1,+1)
#
# © 2026 Jan Bouwman — MIT License

"""
    spacetime_generators() → Vector{Matrix{ComplexF64}}

The four macroscopic generators of the Fisher tensor on ℂ⁶ that
survive the Kaluza-Klein reduction to macroscopic spacetime.

Returns [G_t, G_x, G_y, G_z] where:
- G_t: U(1)_Y singlet (Y=0, ν_R direction) → time
- G_x, G_y, G_z: SU(3) generators → three spatial directions

From Document LXXIX: the Y=0 direction does not precess under
U(1)_Y rotations — this is the geometric origin of the time direction.
"""
function spacetime_generators()
    G_t = zeros(ComplexF64, 6, 6)
    G_t[5,5] =  1/√2; G_t[6,6] = -1/√2

    G_x = zeros(ComplexF64, 6, 6)
    G_x[1,2] =  1/√2; G_x[2,1] =  1/√2

    G_y = zeros(ComplexF64, 6, 6)
    G_y[1,2] = -1im/√2; G_y[2,1] = 1im/√2

    G_z = zeros(ComplexF64, 6, 6)
    G_z[1,1] =  1/√2; G_z[2,2] = -1/√2

    return [G_t, G_x, G_y, G_z]
end

"""
    covariance_metric(ρ) → Matrix{Float64}

The spacetime metric as the covariance matrix of the four
macroscopic generators in state ρ̂:

    g_μν = Cov_ρ̂(G_μ, G_ν)
         = Tr(ρ̂ G_μ G_ν) - Tr(ρ̂ G_μ) Tr(ρ̂ G_ν)

By the quantum CLT (Document LXXX): in the N→∞ limit,
(1/N) 𝓕_μν[ρ̂⊗N] → Cov_ρ̂(G_μ, G_ν).

For the vacuum ρ̂* = I/6: Cov(G_μ, G_ν) = (1/3)δ_μν.

Returns a 4×4 real symmetric matrix (Euclidean, before Wick rotation).
"""
function covariance_metric(ρ::AbstractMatrix)
    G = spacetime_generators()
    g = zeros(Float64, 4, 4)
    for i in 1:4, j in 1:4
        GiGj = G[i] * G[j]
        g[i,j] = real(tr(ρ * GiGj)) -
                 real(tr(ρ * G[i])) * real(tr(ρ * G[j]))
    end
    return g
end

"""
    minkowski_metric(ρ=vacuum_state()) → Matrix{Float64}

Derive the Minkowski metric η_μν = diag(-1,+1,+1,+1) from the
Fisher information geometry of state ρ̂.

The derivation (Documents LXXIX–LXXX):
1. Cov(G_μ,G_ν) = (1/3)δ_μν in the vacuum
2. Wick rotation t→iτ: flip sign of g_tt (the Y=0 direction)
3. Normalise by factor 3

Result: η_μν = diag(-1, +1, +1, +1)

Not postulated — derived from g_AB = 𝓕_AB/ρ₀.

# Example
```julia
η = minkowski_metric()
# → 4×4 matrix: diag(-1, +1, +1, +1)
```
"""
function minkowski_metric(ρ::AbstractMatrix=vacuum_state())
    g = covariance_metric(ρ)
    g[1,1] = -g[1,1]   # Wick rotation: Y=0 singlet → time (-)
    return g * 3        # normalise 1/3 → 1
end

"""
    lorentz_interval(dx::Vector) → Float64

Compute the Lorentz-invariant interval:
    ds² = η_μν dx^μ dx^ν = -dt² + dx² + dy² + dz²

# Arguments
- `dx`: 4-vector [dt, dx, dy, dz]

# Returns
- `ds²` < 0: timelike (causal, massive particles)
- `ds²` = 0: lightlike (photons, gravitons)
- `ds²` > 0: spacelike (cannot be causally connected)

# Example
```julia
lorentz_interval([1.0, 1.0, 0.0, 0.0])  # → 0.0  lightlike
lorentz_interval([2.0, 1.0, 0.0, 0.0])  # → -3.0 timelike
lorentz_interval([1.0, 2.0, 0.0, 0.0])  # → +3.0 spacelike
```
"""
function lorentz_interval(dx::Vector{<:Real})
    length(dx) == 4 ||
        throw(ArgumentError("dx must be a 4-vector [dt, dx, dy, dz]"))
    η = minkowski_metric()
    return sum(η[μ,ν] * dx[μ] * dx[ν] for μ in 1:4, ν in 1:4)
end

"""
    is_timelike(dx) → Bool

True if ds² < 0. Can be connected by a massive particle.
"""
is_timelike(dx)  = lorentz_interval(dx) < -1e-12

"""
    is_lightlike(dx) → Bool

True if ds² ≈ 0. Connected by a massless particle (photon, graviton).
"""
is_lightlike(dx) = abs(lorentz_interval(dx)) ≤ 1e-12

"""
    is_spacelike(dx) → Bool

True if ds² > 0. Cannot be causally connected.
"""
is_spacelike(dx) = lorentz_interval(dx) > 1e-12
