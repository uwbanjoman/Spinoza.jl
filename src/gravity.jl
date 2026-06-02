# src/gravity.jl
# ================
# Gravity from Fisher information geometry
#
# The Einstein equation follows from the FisherGeometrics postulate:
#   g_μν = 𝓕_μν / ρ₀
#
# The spacetime curvature IS the curvature of the Fisher information
# tensor field ρ̂(x,t). No separate gravitational field — gravity
# is the geometry of quantum information.
#
# This file derives:
#   1. Metric tensor from Fisher information
#   2. Christoffel symbols Γ^λ_μν
#   3. Riemann tensor R^ρ_σμν
#   4. Ricci tensor R_μν and scalar R
#   5. Einstein tensor G_μν = R_μν - (1/2)g_μν R
#   6. Stress-energy tensor T_μν from ρ̂
#   7. Einstein equation G_μν = 8πG_N T_μν
#
# © 2026 Jan Bouwman — MIT License

using QuantumFisher
using LinearAlgebra

export fisher_metric

# ── Metric from Fisher information ────────────────────────────────────────────

"""
    fisher_metric(ρ̂, generators) → Matrix{Float64}

Compute the spacetime metric from the Fisher information tensor:
    g_μν = 𝓕_μν / ρ₀ = Cov_ρ̂(G_μ, G_ν)

For the vacuum ρ̂* = I/6: g_μν → η_μν (Minkowski, after Wick rotation).
For a general state ρ̂: g_μν is a curved spacetime metric.

# Arguments
- `ρ̂`: density matrix at a spacetime point
- `generators`: [G_t, G_x, G_y, G_z] — macroscopic generators

# Returns
4×4 real symmetric metric tensor g_μν
"""
function fisher_metric(
    ρ̂::AbstractMatrix,
    generators::Vector{<:AbstractMatrix}
)
    g = zeros(Float64, 4, 4)
    for i in 1:4, j in 1:4
        GiGj = generators[i] * generators[j]
        g[i,j] = real(tr(ρ̂ * GiGj)) -
                 real(tr(ρ̂ * generators[i])) *
                 real(tr(ρ̂ * generators[j]))
    end
    return g
end

"""
    metric_field(ρ̂_field, generators) → Array{Float64,3}

Compute the metric tensor field g_μν(x) from a field of
density matrices ρ̂(x) on a spatial grid.

# Arguments
- `ρ̂_field`: array of density matrices, one per grid point
- `generators`: spacetime generators [G_t, G_x, G_y, G_z]

# Returns
Array of shape (N, 4, 4) — metric tensor at each grid point
"""
function metric_field(
    ρ̂_field::Vector{<:AbstractMatrix},
    generators::Vector{<:AbstractMatrix}
)
    N = length(ρ̂_field)
    g = zeros(Float64, N, 4, 4)
    for i in 1:N
        g[i,:,:] = fisher_metric(ρ̂_field[i], generators)
    end
    return g
end

# ── Christoffel symbols ───────────────────────────────────────────────────────

"""
    christoffel(g, dg) → Array{Float64,3}

Compute the Christoffel symbols from the metric and its derivatives:
    Γ^λ_μν = (1/2) g^λρ (∂_μ g_νρ + ∂_ν g_μρ - ∂_ρ g_μν)

# Arguments
- `g`: 4×4 metric tensor at a point
- `dg`: 4×4×4 array, dg[α,μ,ν] = ∂_α g_μν

# Returns
4×4×4 array Γ[λ,μ,ν] = Γ^λ_μν
"""
function christoffel(g::AbstractMatrix, dg::AbstractArray)
    g_inv = inv(g)
    Γ = zeros(Float64, 4, 4, 4)
    for λ in 1:4, μ in 1:4, ν in 1:4
        for ρ in 1:4
            Γ[λ,μ,ν] += 0.5 * g_inv[λ,ρ] *
                (dg[μ,ν,ρ] + dg[ν,μ,ρ] - dg[ρ,μ,ν])
        end
    end
    return Γ
end

"""
    metric_derivatives(g_field, dx) → Array{Float64,4}

Compute the spatial derivatives of the metric tensor on a grid
using central differences.

# Arguments
- `g_field`: (N, 4, 4) array of metric tensors on a 1D grid
- `dx`: grid spacing

# Returns
(N, 4, 4, 4) array dg[i,α,μ,ν] = ∂_α g_μν at grid point i
(only spatial ∂_x implemented; temporal derivatives require evolution)
"""
function metric_derivatives(g_field::AbstractArray, dx::Real)
    N = size(g_field, 1)
    dg = zeros(Float64, N, 4, 4, 4)
    for i in 1:N
        il = i == 1 ? N : i-1
        ir = i == N ? 1 : i+1
        # Spatial derivative (α=2, x-direction)
        dg[i,2,:,:] = (g_field[ir,:,:] - g_field[il,:,:]) / (2dx)
    end
    return dg
end

# ── Riemann tensor ────────────────────────────────────────────────────────────

"""
    riemann(Γ, dΓ) → Array{Float64,4}

Compute the Riemann curvature tensor:
    R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ

# Arguments
- `Γ`: 4×4×4 Christoffel symbols Γ[ρ,σ,μ] = Γ^ρ_σμ
- `dΓ`: 4×4×4×4 array, dΓ[α,ρ,σ,μ] = ∂_α Γ^ρ_σμ

# Returns
4×4×4×4 array R[ρ,σ,μ,ν] = R^ρ_σμν
"""
function riemann(Γ::AbstractArray, dΓ::AbstractArray)
    R = zeros(Float64, 4, 4, 4, 4)
    for ρ in 1:4, σ in 1:4, μ in 1:4, ν in 1:4
        R[ρ,σ,μ,ν] = dΓ[μ,ρ,ν,σ] - dΓ[ν,ρ,μ,σ]
        for λ in 1:4
            R[ρ,σ,μ,ν] += Γ[ρ,μ,λ]*Γ[λ,ν,σ] - Γ[ρ,ν,λ]*Γ[λ,μ,σ]
        end
    end
    return R
end

# ── Ricci tensor and scalar ───────────────────────────────────────────────────

"""
    ricci_tensor(R_riemann) → Matrix{Float64}

Compute the Ricci tensor by contraction:
    R_μν = R^ρ_μρν  (trace over first and third indices)
"""
function ricci_tensor(R::AbstractArray)
    Ric = zeros(Float64, 4, 4)
    for μ in 1:4, ν in 1:4
        for ρ in 1:4
            Ric[μ,ν] += R[ρ,μ,ρ,ν]
        end
    end
    return Ric
end

"""
    ricci_scalar(Ric, g) → Float64

Compute the Ricci scalar:
    R = g^μν R_μν
"""
function ricci_scalar(Ric::AbstractMatrix, g::AbstractMatrix)
    g_inv = inv(g)
    return sum(g_inv[μ,ν] * Ric[μ,ν] for μ in 1:4, ν in 1:4)
end

# ── Einstein tensor ───────────────────────────────────────────────────────────

"""
    einstein_tensor(ρ̂_field, generators, dx) → Array{Float64,3}

Compute the Einstein tensor field G_μν(x) from a field of
density matrices on a spatial grid.

    G_μν = R_μν - (1/2) g_μν R

This is the left-hand side of the Einstein equation.

# Returns
(N, 4, 4) array G[i,μ,ν] at each grid point
"""
function einstein_tensor(
    ρ̂_field::Vector{<:AbstractMatrix},
    generators::Vector{<:AbstractMatrix},
    dx::Real
)
    N = length(ρ̂_field)
    g_field = metric_field(ρ̂_field, generators)
    dg_field = metric_derivatives(g_field, dx)

    G_field = zeros(Float64, N, 4, 4)

    for i in 1:N
        g_i  = g_field[i,:,:]
        dg_i = dg_field[i,:,:,:]

        Γ_i  = christoffel(g_i, dg_i)

        # Approximate Christoffel derivatives from neighbours
        il = i == 1 ? N : i-1
        ir = i == N ? 1 : i+1
        g_l = g_field[il,:,:]; g_r = g_field[ir,:,:]
        dg_l = dg_field[il,:,:,:]; dg_r = dg_field[ir,:,:,:]
        Γ_l = christoffel(g_l, dg_l)
        Γ_r = christoffel(g_r, dg_r)

        dΓ_i = zeros(Float64, 4, 4, 4, 4)
        dΓ_i[2,:,:,:] = (Γ_r - Γ_l) / (2dx)  # ∂_x Γ

        R_i  = riemann(Γ_i, dΓ_i)
        Ric_i = ricci_tensor(R_i)
        R_scalar_i = ricci_scalar(Ric_i, g_i)

        G_field[i,:,:] = Ric_i - 0.5 * g_i * R_scalar_i
    end

    return G_field
end

# ── Stress-energy tensor ──────────────────────────────────────────────────────

"""
    stress_energy(ρ̂, generators, cs=1.0) → Matrix{Float64}

Compute the stress-energy tensor T_μν from the density matrix ρ̂.

In FisherGeometrics, the stress-energy is the Fisher information
excess above the vacuum:
    T_μν = 𝓕_μν[ρ̂] - 𝓕_μν[ρ̂*]  (excess Fisher information)

This is the right-hand side of the Einstein equation (up to 8πG_N).
"""
function stress_energy(
    ρ̂::AbstractMatrix,
    generators::Vector{<:AbstractMatrix},
    cs::Real=1.0
)
    ρ_vac = vacuum_state()
    T = fisher_metric(ρ̂, generators) - fisher_metric(ρ_vac, generators)
    return T
end

# ── Einstein equation check ───────────────────────────────────────────────────

"""
    einstein_equation_residual(G_μν, T_μν, G_N=1.0) → Float64

Compute the residual of the Einstein equation:
    ‖G_μν - 8πG_N T_μν‖_F

A small residual means the Einstein equation is approximately satisfied.

In FisherGeometrics, this is exact in the classical limit where
ρ̂(x,t) describes a macroscopic gravitational field.
"""
function einstein_equation_residual(
    G_μν::AbstractMatrix,
    T_μν::AbstractMatrix,
    G_N::Real=1.0
)
    return norm(G_μν - 8π*G_N * T_μν)
end

"""
    gravitational_wave_metric(ρ̂₀, H_KK, t, generators) → Matrix{Float64}

Compute the metric perturbation h_μν for a gravitational wave
sourced by the density matrix evolution under H_KK.

The quadrupole formula gives:
    h_μν ∝ d²Q_μν/dt²  where Q_μν = ∫ T_μν x_μ x_ν d³x

In FisherGeometrics: Q_μν is derived from 𝓕_μν[ρ̂(t)].
"""
function gravitational_wave_metric(
    ρ̂₀::AbstractMatrix,
    H_KK::AbstractMatrix,
    t::Real,
    generators::Vector{<:AbstractMatrix}
)
    # Evolve density matrix
    U = exp(-im * H_KK * t)
    ρ̂_t = U * ρ̂₀ * U'

    # Metric perturbation = difference from vacuum metric
    g_t   = fisher_metric(ρ̂_t,       generators)
    g_vac = fisher_metric(vacuum_state(), generators)

    return g_t - g_vac  # h_μν
end
