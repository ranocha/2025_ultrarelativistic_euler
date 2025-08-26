# Install packages
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Load packages
using LinearAlgebra: normalize, norm
using Test: @test, @testset

using Trixi
using Trixi: sqrt # return NaN for negative values instead of throwing an error
using Trixi: density # GLMakie exports density, too
using OrdinaryDiffEqCore: PIDController
using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK

using HDF5


const DATA_DIR = joinpath(dirname(@__DIR__), "data")
if !isdir(DATA_DIR)
    mkdir(DATA_DIR)
end


#####################################################################
# Implementation of the ultra-relativistic Euler equations in 2D
struct UltraRelativisticEulerEquations2D <: Trixi.AbstractEquations{2, 3}
end

function Trixi.varnames(::typeof(cons2cons),  ::UltraRelativisticEulerEquations2D)
    ("rho_v1", "rho_v2", "rho_e")
end

function Trixi.varnames(::typeof(cons2prim),  ::UltraRelativisticEulerEquations2D)
    ("v1", "v2", "p")
end

# Convert conservative variables to primitive variables
@inline function Trixi.cons2prim(u,
                                 equations::UltraRelativisticEulerEquations2D)
    rho_v1, rho_v2, rho_e = u

    p = (sqrt(4 * rho_e^2 - 3 * (rho_v1^2 + rho_v2^2)) - rho_e) / 3
    factor = 1 / sqrt(4 * p * (rho_e + p))
    v1 = rho_v1 * factor
    v2 = rho_v2 * factor

    return SVector(v1, v2, p)
end

# Convert primitive variables to conservative variables
@inline function Trixi.prim2cons(prim,
                                 equations::UltraRelativisticEulerEquations2D)
    v1, v2, p = prim

    factor = 4 * p * sqrt(1 + v1^2 + v2^2)
    rho_v1 = v1 * factor
    rho_v2 = v2 * factor
    rho_e = p * (3 + 4 * (v1^2 + v2^2))

    return SVector(rho_v1, rho_v2, rho_e)
end

@inline function Trixi.pressure(u,
                                equations::UltraRelativisticEulerEquations2D)
    rho_v1, rho_v2, rho_e = u

    p = (sqrt(4 * rho_e^2 - 3 * (rho_v1^2 + rho_v2^2)) - rho_e) / 3

    return p
end

@inline function pressure_velocity(u, equations::UltraRelativisticEulerEquations2D)
    v1, v2, p = cons2prim(u, equations)
    return p * sqrt(1 + v1^2 + v2^2)
end

# Calculate normal flux for a single point
@inline function Trixi.flux(u, orientation::Integer,
                            equations::UltraRelativisticEulerEquations2D)
    v1, v2, p = cons2prim(u, equations)
    factor = sqrt(1 + v1^2 + v2^2)

    if orientation == 1 # x-direction
        f1 = 4 * p * v1^2 + p
        f2 = 4 * p * v1 * v2
        fe = 4 * p * v1 * factor
        return SVector(f1, f2, fe)
    else # y-direction
        f1 = 4 * p * v1 * v2
        f2 = 4 * p * v2^2 + p
        fe = 4 * p * v2 * factor
        return SVector(f1, f2, fe)
    end
end

@inline function Trixi.flux(u, normal_direction,
                            equations::UltraRelativisticEulerEquations2D)
    v1, v2, p = cons2prim(u, equations)
    factor = sqrt(1 + v1^2 + v2^2)

    v_dot_n = (v1 * normal_direction[1] +
               v2 * normal_direction[2])

    f1 = 4 * p * v1 * v_dot_n + p * normal_direction[1]
    f2 = 4 * p * v2 * v_dot_n + p * normal_direction[2]
    fe = 4 * p * v_dot_n * factor

    return SVector(f1, f2, fe)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr,
                                           orientation::Integer,
                                           equations::UltraRelativisticEulerEquations2D)
    v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    else # orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    end

    denom_ll = 3 + 2 * (v1_ll^2 + v2_ll^2)
    denom_rr = 3 + 2 * (v1_rr^2 + v2_rr^2)

    lambda_ll = (2 * v_ll^2 + sqrt(2 * (v1_ll^2 + v2_ll^2 - v_ll^2) + 3)) / denom_ll
    lambda_rr = (2 * v_rr^2 + sqrt(2 * (v1_rr^2 + v2_rr^2 - v_rr^2) + 3)) / denom_rr

    λ_max = max(lambda_ll, lambda_rr)
    return λ_max
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr,
                                           normal_direction,
                                           equations::UltraRelativisticEulerEquations2D)
    v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    norm_ = norm(normal_direction)

    v_dot_n_ll = (v1_ll * normal_direction[1] +
                  v2_ll * normal_direction[2]) / norm_
    v_dot_n_rr = (v1_rr * normal_direction[1] +
                  v2_rr * normal_direction[2]) / norm_

    denom_ll = 3 + 2 * (v1_ll^2 + v2_ll^2)
    denom_rr = 3 + 2 * (v1_rr^2 + v2_rr^2)

    lambda_ll = (2 * v_dot_n_ll^2 + sqrt(2 * (v1_ll^2 + v2_ll^2 - v_dot_n_ll^2) + 3)) / denom_ll
    lambda_rr = (2 * v_dot_n_rr^2 + sqrt(2 * (v1_rr^2 + v2_rr^2 - v_dot_n_rr^2) + 3)) / denom_rr

    λ_max = max(lambda_ll, lambda_rr) * norm_
    return λ_max
end

@inline function Trixi.max_abs_speeds(u, equations::UltraRelativisticEulerEquations2D)
    v1, v2, p = cons2prim(u, equations)

    denom = 3 + 2 * (v1^2 + v2^2)

    lambda_1 = (2 * v1^2 + sqrt(2 * v2^2 + 3)) / denom
    lambda_2 = (2 * v2^2 + sqrt(2 * v1^2 + 3)) / denom

    return lambda_1, lambda_2
end

# Convert conservative variables to entropy variables
@inline function Trixi.cons2entropy(u,
                                    equations::UltraRelativisticEulerEquations2D)
    v1, v2, p = cons2prim(u, equations)

    # compute the fourth root ("tessaract root") of the pressure
    sqrt_p = sqrt(p)
    tert_p = sqrt(sqrt_p)

    w1 = -v1 / (4 * tert_p)
    w2 = -v2 / (4 * tert_p)
    we = sqrt(1 + v1^2 + v2^2) / (4 * tert_p)

    return SVector(w1, w2, we)
end

@inline function flux_thein_ranocha(u_ll, u_rr,
                                    orientation::Integer,
                                    equations::UltraRelativisticEulerEquations2D)
    v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # compute the fourth roots ("tessaract roots") of the pressure
    sqrt_p_ll = sqrt(p_ll)
    tert_p_ll = sqrt(sqrt_p_ll)
    sqrt_p_rr = sqrt(p_rr)
    tert_p_rr = sqrt(sqrt_p_rr)

    if orientation == 1
        f1 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr)^2 + (p_ll + p_rr))
        f2 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr))
        fe = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (sqrt(1 + v1_ll^2 + v2_ll^2) / tert_p_ll + sqrt(1 + v1_rr^2 + v2_rr^2) / tert_p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr))
    else # orientation == 2
        f1 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr))
        f2 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr)^2 + (p_ll + p_rr))
        fe = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (sqrt(1 + v1_ll^2 + v2_ll^2) / tert_p_ll + sqrt(1 + v1_rr^2 + v2_rr^2) / tert_p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr))
    end

    return SVector(f1, f2, fe)
end

@inline function flux_thein_ranocha(u_ll, u_rr,
                                    normal_direction,
                                    equations::UltraRelativisticEulerEquations2D)
    v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # compute the fourth roots ("tessaract roots") of the pressure
    sqrt_p_ll = sqrt(p_ll)
    tert_p_ll = sqrt(sqrt_p_ll)
    sqrt_p_rr = sqrt(p_rr)
    tert_p_rr = sqrt(sqrt_p_rr)

    v_dot_n_ll = (v1_ll * normal_direction[1] +
                  v2_ll * normal_direction[2])
    v_dot_n_rr = (v1_rr * normal_direction[1] +
                  v2_rr * normal_direction[2])

    f1 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                  (v1_ll / tert_p_ll + v1_rr / tert_p_rr) *
                  (v_dot_n_ll / tert_p_ll + v_dot_n_rr / tert_p_rr)
                  + (p_ll + p_rr) * normal_direction[1])
    f2 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                  (v2_ll / tert_p_ll + v2_rr / tert_p_rr) *
                  (v_dot_n_ll / tert_p_ll + v_dot_n_rr / tert_p_rr)
                  + (p_ll + p_rr) * normal_direction[2])
    fe = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                  (sqrt(1 + v1_ll^2 + v2_ll^2) / tert_p_ll + sqrt(1 + v1_rr^2 + v2_rr^2) / tert_p_rr) *
                  (v_dot_n_ll / tert_p_ll + v_dot_n_rr / tert_p_rr))

    return SVector(f1, f2, fe)
end

function unit_tests_2d()
@testset "UltraRelativisticEulerEquations2D" begin
    equations = UltraRelativisticEulerEquations2D()

    u_ll = SVector(5.5, -0.5, 7.0)
    u_rr = SVector(-3.0, 2.0, 8.0)

    orientations = 1:2
    normal_directions = [SVector(1.0, 0.0),
                         SVector(0.0, 1.0),
                         SVector(1.0, 1.0),
                         SVector(1.0, -1.0),
                         SVector(-1.0, 1.0),
                         SVector(-1.0, -1.0)]
    factors = (0.12345, π / 3, -exp(1), 1.23456789)
    fluxes = (flux_thein_ranocha, flux_lax_friedrichs)

    @testset "Consistency checks of numerical fluxes" begin
        @testset "$fnum" for fnum in fluxes
            @testset "orientation $orientation" for orientation in orientations
                f = fnum(u_ll, u_ll, orientation, equations)
                @test f ≈ flux(u_ll, orientation, equations)

                f = fnum(u_rr, u_rr, orientation, equations)
                @test f ≈ flux(u_rr, orientation, equations)
            end
        end
        @testset "$fnum" for fnum in fluxes
            @testset "normal_direction $normal_direction" for normal_direction in normal_directions
                f = fnum(u_ll, u_ll, normal_direction, equations)
                @test f ≈ flux(u_ll, normal_direction, equations)

                f = fnum(u_rr, u_rr, normal_direction, equations)
                @test f ≈ flux(u_rr, normal_direction, equations)
            end
        end
    end

    @testset "Scaling of numerical fluxes by normal_direction" begin
        @testset "factor $factor" for factor in factors
            for normal_direction in normal_directions
                f = flux_thein_ranocha(u_ll, u_rr, normal_direction, equations)
                f_scaled = flux_thein_ranocha(u_ll, u_rr, factor * normal_direction, equations)
                @test f_scaled ≈ factor * f
            end
        end
    end
end
end


#####################################################################
# Implementation of the ultra-relativistic Euler equations in 3D
struct UltraRelativisticEulerEquations3D <: Trixi.AbstractEquations{3, 4}
end

function Trixi.varnames(::typeof(cons2cons),  ::UltraRelativisticEulerEquations3D)
    ("rho_v1", "rho_v2", "rho_v3", "rho_e")
end

function Trixi.varnames(::typeof(cons2prim),  ::UltraRelativisticEulerEquations3D)
    ("v1", "v2", "v3", "p")
end

# Convert conservative variables to primitive variables
@inline function Trixi.cons2prim(u,
                                 equations::UltraRelativisticEulerEquations3D)
    rho_v1, rho_v2, rho_v3, rho_e = u

    p = (sqrt(4 * rho_e^2 - 3 * (rho_v1^2 + rho_v2^2 + rho_v3^2)) - rho_e) / 3
    factor = 1 / sqrt(4 * p * (rho_e + p))
    v1 = rho_v1 * factor
    v2 = rho_v2 * factor
    v3 = rho_v3 * factor

    return SVector(v1, v2, v3, p)
end

# Convert primitive variables to conservative variables
@inline function Trixi.prim2cons(prim,
                                 equations::UltraRelativisticEulerEquations3D)
    v1, v2, v3, p = prim

    factor = 4 * p * sqrt(1 + v1^2 + v2^2 + v3^2)
    rho_v1 = v1 * factor
    rho_v2 = v2 * factor
    rho_v3 = v3 * factor
    rho_e = p * (3 + 4 * (v1^2 + v2^2 + v3^2))

    return SVector(rho_v1, rho_v2, rho_v3, rho_e)
end

@inline function Trixi.pressure(u,
                                equations::UltraRelativisticEulerEquations3D)
    rho_v1, rho_v2, rho_v3, rho_e = u

    p = (sqrt(4 * rho_e^2 - 3 * (rho_v1^2 + rho_v2^2 + rho_v3^2)) - rho_e) / 3

    return p
end

@inline function pressure_velocity(u, equations::UltraRelativisticEulerEquations3D)
    v1, v2, v3, p = cons2prim(u, equations)
    return p * sqrt(1 + v1^2 + v2^2 + v3^2)
end

# Calculate normal flux for a single point
@inline function Trixi.flux(u, orientation::Integer,
                            equations::UltraRelativisticEulerEquations3D)
    v1, v2, v3, p = cons2prim(u, equations)
    factor = sqrt(1 + v1^2 + v2^2 + v3^2)

    if orientation == 1 # x-direction
        f1 = 4 * p * v1^2 + p
        f2 = 4 * p * v1 * v2
        f3 = 4 * p * v1 * v3
        fe = 4 * p * v1 * factor
        return SVector(f1, f2, f3, fe)
    elseif orientation == 2 # y-direction
        f1 = 4 * p * v1 * v2
        f2 = 4 * p * v2^2 + p
        f3 = 4 * p * v2 * v3
        fe = 4 * p * v2 * factor
        return SVector(f1, f2, f3, fe)
    else # z-direction
        f1 = 4 * p * v1 * v3
        f2 = 4 * p * v2 * v3
        f3 = 4 * p * v3^2 + p
        fe = 4 * p * v3 * factor
        return SVector(f1, f2, f3, fe)
    end
end

@inline function Trixi.flux(u, normal_direction,
                            equations::UltraRelativisticEulerEquations3D)
    v1, v2, v3, p = cons2prim(u, equations)
    factor = sqrt(1 + v1^2 + v2^2 + v3^2)
    v_dot_n = (v1 * normal_direction[1] +
               v2 * normal_direction[2] +
               v3 * normal_direction[3])

    f1 = 4 * p * v1 * v_dot_n + p * normal_direction[1]
    f2 = 4 * p * v2 * v_dot_n + p * normal_direction[2]
    f3 = 4 * p * v3 * v_dot_n + p * normal_direction[3]
    fe = 4 * p * v_dot_n * factor

    return SVector(f1, f2, f3, fe)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function Trixi.max_abs_speed_naive(u_ll, u_rr,
                                           orientation::Integer,
                                           equations::UltraRelativisticEulerEquations3D)
    v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    elseif orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    else # orientation == 3
        v_ll = v3_ll
        v_rr = v3_rr
    end

    denom_ll = 3 + 2 * (v1_ll^2 + v2_ll^2 + v3_ll^2)
    denom_rr = 3 + 2 * (v1_rr^2 + v2_rr^2 + v3_rr^2)

    lambda_ll = (2 * v_ll^2 + sqrt(2 * (v1_ll^2 + v2_ll^2 + v3_ll^2 - v_ll^2) + 3)) / denom_ll
    lambda_rr = (2 * v_rr^2 + sqrt(2 * (v1_rr^2 + v2_rr^2 + v3_rr^2 - v_rr^2) + 3)) / denom_rr

    λ_max = max(lambda_ll, lambda_rr)
    return λ_max
end

@inline function Trixi.max_abs_speed_naive(u_ll, u_rr,
                                           normal_direction,
                                           equations::UltraRelativisticEulerEquations3D)
    v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    norm_ = norm(normal_direction)

    v_dot_n_ll = (v1_ll * normal_direction[1] +
                  v2_ll * normal_direction[2] +
                  v3_ll * normal_direction[3]) / norm_
    v_dot_n_rr = (v1_rr * normal_direction[1] +
                  v2_rr * normal_direction[2] +
                  v3_rr * normal_direction[3]) / norm_

    denom_ll = 3 + 2 * (v1_ll^2 + v2_ll^2 + v3_ll^2)
    denom_rr = 3 + 2 * (v1_rr^2 + v2_rr^2 + v3_rr^2)

    lambda_ll = (2 * v_dot_n_ll^2 + sqrt(2 * (v1_ll^2 + v2_ll^2 + v3_ll^2 - v_dot_n_ll^2) + 3)) / denom_ll
    lambda_rr = (2 * v_dot_n_rr^2 + sqrt(2 * (v1_rr^2 + v2_rr^2 + v3_rr^2 - v_dot_n_rr^2) + 3)) / denom_rr

    λ_max = max(lambda_ll, lambda_rr) * norm_
    return λ_max
end

@inline function Trixi.max_abs_speeds(u, equations::UltraRelativisticEulerEquations3D)
    v1, v2, v3, p = cons2prim(u, equations)

    denom = 3 + 2 * (v1^2 + v2^2 + v3^2)

    lambda_1 = (2 * v1^2 + sqrt(2 * (v2^2 + v3^2) + 3)) / denom
    lambda_2 = (2 * v2^2 + sqrt(2 * (v1^2 + v3^2) + 3)) / denom
    lambda_3 = (2 * v3^2 + sqrt(2 * (v1^2 + v2^2) + 3)) / denom

    return lambda_1, lambda_2, lambda_3
end

# Convert conservative variables to entropy variables
@inline function Trixi.cons2entropy(u,
                                    equations::UltraRelativisticEulerEquations3D)
    v1, v2, v3, p = cons2prim(u, equations)

    # compute the fourth root ("tessaract root") of the pressure
    sqrt_p = sqrt(p)
    tert_p = sqrt(sqrt_p)

    w1 = -v1 / (4 * tert_p)
    w2 = -v2 / (4 * tert_p)
    w3 = -v3 / (4 * tert_p)
    we = sqrt(1 + v1^2 + v2^2 + v3^2) / (4 * tert_p)

    return SVector(w1, w2, w3, we)
end

@inline function flux_thein_ranocha(u_ll, u_rr,
                                    orientation::Integer,
                                    equations::UltraRelativisticEulerEquations3D)
    v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    # compute the fourth roots ("tessaract roots") of the pressure
    sqrt_p_ll = sqrt(p_ll)
    tert_p_ll = sqrt(sqrt_p_ll)
    sqrt_p_rr = sqrt(p_rr)
    tert_p_rr = sqrt(sqrt_p_rr)

    if orientation == 1
        f1 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr)^2 + (p_ll + p_rr))
        f2 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr))
        f3 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr) *
                      (v3_ll / tert_p_ll + v3_rr / tert_p_rr))
        fe = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (sqrt(1 + v1_ll^2 + v2_ll^2 + v3_ll^2) / tert_p_ll + sqrt(1 + v1_rr^2 + v2_rr^2 + v3_rr^2) / tert_p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr))
    elseif orientation == 2
        f1 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr))
        f2 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr)^2 + (p_ll + p_rr))
        f3 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr) *
                      (v3_ll / tert_p_ll + v3_rr / tert_p_rr))
        fe = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (sqrt(1 + v1_ll^2 + v2_ll^2 + v3_ll^2) / tert_p_ll + sqrt(1 + v1_rr^2 + v2_rr^2 + v3_rr^2) / tert_p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr))
    else # orientation == 3
        f1 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v1_ll / tert_p_ll + v1_rr / tert_p_rr) *
                      (v3_ll / tert_p_ll + v3_rr / tert_p_rr))
        f2 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v2_ll / tert_p_ll + v2_rr / tert_p_rr) *
                      (v3_ll / tert_p_ll + v3_rr / tert_p_rr))
        f3 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (v3_ll / tert_p_ll + v3_rr / tert_p_rr)^2 + (p_ll + p_rr))
        fe = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                      (sqrt(1 + v1_ll^2 + v2_ll^2 + v3_ll^2) / tert_p_ll + sqrt(1 + v1_rr^2 + v2_rr^2 + v3_rr^2) / tert_p_rr) *
                      (v3_ll / tert_p_ll + v3_rr / tert_p_rr))
    end

    return SVector(f1, f2, f3, fe)
end

@inline function flux_thein_ranocha(u_ll, u_rr,
                                    normal_direction,
                                    equations::UltraRelativisticEulerEquations3D)
    v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    # compute the fourth roots ("tessaract roots") of the pressure
    sqrt_p_ll = sqrt(p_ll)
    tert_p_ll = sqrt(sqrt_p_ll)
    sqrt_p_rr = sqrt(p_rr)
    tert_p_rr = sqrt(sqrt_p_rr)

    v_dot_n_ll = (v1_ll * normal_direction[1] +
                  v2_ll * normal_direction[2] +
                  v3_ll * normal_direction[3])
    v_dot_n_rr = (v1_rr * normal_direction[1] +
                  v2_rr * normal_direction[2] +
                  v3_rr * normal_direction[3])

    f1 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                  (v1_ll / tert_p_ll + v1_rr / tert_p_rr) *
                  (v_dot_n_ll / tert_p_ll + v_dot_n_rr / tert_p_rr)
                  + (p_ll + p_rr) * normal_direction[1])
    f2 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                  (v2_ll / tert_p_ll + v2_rr / tert_p_rr) *
                  (v_dot_n_ll / tert_p_ll + v_dot_n_rr / tert_p_rr)
                  + (p_ll + p_rr) * normal_direction[2])
    f3 = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                  (v3_ll / tert_p_ll + v3_rr / tert_p_rr) *
                  (v_dot_n_ll / tert_p_ll + v_dot_n_rr / tert_p_rr)
                  + (p_ll + p_rr) * normal_direction[3])
    fe = 0.5f0 * ((p_ll * sqrt_p_rr + sqrt_p_ll * p_rr) *
                  (sqrt(1 + v1_ll^2 + v2_ll^2 + v3_ll^2) / tert_p_ll + sqrt(1 + v1_rr^2 + v2_rr^2 + v3_rr^2) / tert_p_rr) *
                  (v_dot_n_ll / tert_p_ll + v_dot_n_rr / tert_p_rr))

    return SVector(f1, f2, f3, fe)
end

function unit_tests_3d()
@testset "UltraRelativisticEulerEquations3D" begin
    equations = UltraRelativisticEulerEquations3D()

    u_ll = SVector(5.5, -0.5, -1.5, 7.0)
    u_rr = SVector(-3.0, 2.0, 1.0, 8.0)

    orientations = 1:3
    normal_directions = [SVector(1.0, 0.0, 0.0),
                         SVector(0.0, 1.0, 0.0),
                         SVector(0.0, 0.0, 1.0),
                         SVector(1.0, -1.0, 0.0),
                         SVector(1.0, 0.0, -1.0),
                         SVector(0.0, -1.0, 1.0),
                         SVector(1.0, 1.0, 1.0),
                         SVector(-1.0, 1.0, 1.0),
                         SVector(1.0, -1.0, 1.0),
                         SVector(1.0, 1.0, -1.0)]
    factors = (0.12345, π / 3, -exp(1), 1.23456789)
    fluxes = (flux_thein_ranocha, flux_lax_friedrichs)

    @testset "Consistency checks of numerical fluxes" begin
        @testset "$fnum" for fnum in fluxes
            @testset "orientation $orientation" for orientation in orientations
                f = fnum(u_ll, u_ll, orientation, equations)
                @test f ≈ flux(u_ll, orientation, equations)

                f = fnum(u_rr, u_rr, orientation, equations)
                @test f ≈ flux(u_rr, orientation, equations)
            end
        end
        @testset "$fnum" for fnum in fluxes
            @testset "normal_direction $normal_direction" for normal_direction in normal_directions
                f = fnum(u_ll, u_ll, normal_direction, equations)
                @test f ≈ flux(u_ll, normal_direction, equations)

                f = fnum(u_rr, u_rr, normal_direction, equations)
                @test f ≈ flux(u_rr, normal_direction, equations)
            end
        end
    end

    @testset "Scaling of numerical fluxes by normal_direction" begin
        @testset "factor $factor" for factor in factors
            for normal_direction in normal_directions
                f = flux_thein_ranocha(u_ll, u_rr, normal_direction, equations)
                f_scaled = flux_thein_ranocha(u_ll, u_rr, factor * normal_direction, equations)
                @test f_scaled ≈ factor * f
            end
        end
    end

    @testset "Consistency of Cartesian and curved fluxes" begin
        @testset "$fnum" for fnum in fluxes
            for (orientation, normal_direction) in ((1, SVector(1, 0, 0)),
                                                    (2, SVector(0, 1, 0)),
                                                    (3, SVector(0, 0, 1)))
                f = fnum(u_ll, u_rr, orientation, equations)
                f_curved = fnum(u_ll, u_rr, normal_direction, equations)
                @test f ≈ f_curved
            end
        end
    end
end
end


###########################################################
# Example 1 of Kunik, Kolb, Müller, and Thein (2024)
function save_shock_and_stationary_2d(; base_level = 3,
                                        med_level = 5,
                                        max_level = 9,
                                        amr_interval = 5,
                                        volume_flux = flux_thein_ranocha,
                                        MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations2D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations2D)
        p = 1.0
        if iszero(x)
            v = zero(x)
        else
            v = -normalize(x)
        end
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else # P4estMesh, T8codeMesh
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc)
    end

    coordinates_min = (-2.0, -2.0)
    coordinates_max = (+2.0, +2.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^5,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 0.5,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 1.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 100)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 2, length = 1000) / sqrt(2)
    curve = vcat(x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    data_v = @. (v1 + v2) / sqrt(2)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            append!(data_v, @. (v1 + v2) / sqrt(2))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "shock_and_stationary_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end


# Same setup as Example 1 of Kunik, Kolb, Müller, and Thein (2024),
# but in 3D instead of 2D
function save_shock_and_stationary_3d(; base_level = 3,
                                        med_level = 5,
                                        max_level = 7,
                                        amr_interval = 5,
                                        volume_flux = flux_thein_ranocha,
                                        MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations3D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations3D)
        p = 1.0
        if iszero(x)
            v = zero(x)
        else
            v = -normalize(x)
        end
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else # P4estMesh, T8codeMesh
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc,
                                   :z_neg => bc, :z_pos => bc)
    end

    coordinates_min = (-2.0, -2.0, -2.0)
    coordinates_max = (+2.0, +2.0, +2.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^6,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 0.5,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 1.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 25)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 2, length = 1000) / sqrt(3)
    curve = vcat(x', x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    v3 = view(pd.data, :, 3)
    data_v = @. (v1 + v2 + v3) / sqrt(3)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            v3 = view(pd.data, :, 3)
            append!(data_v, @. (v1 + v2 + v3) / sqrt(3))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "shock_and_stationary_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end


# Similar setup as Example 1 of Kunik, Kolb, Müller, and Thein (2024),
# but with smaller absolute value of the velocity to check entropy
# conservation
function save_check_ec_2d(; initial_refinement_level = 6,
                            MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations2D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations2D)
        p = 1.0
        if iszero(x)
            v = zero(x)
        else
            v = -0.2 * normalize(x)
        end
        return prim2cons(SVector(v..., p), equations)
    end

    coordinates_min = (-2.0, -2.0)
    coordinates_max = (+2.0, +2.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level,
                        n_cells_max = 10^5,
                        periodicity = true)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level,
                        periodicity = true)
    end

    volume_flux = flux_thein_ranocha
    solver = DGSEM(polydeg = 3, surface_flux = volume_flux,
                   volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

    tspan = (0.0, 1.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 50)
    callbacks = CallbackSet(summary_callback, alive_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    data_t = Vector{Float64}()
    data_entropy_rate = Vector{Float64}()
    du_ode = similar(integrator.u)

    for _ in integrator
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            Trixi.rhs!(du_ode, integrator.u, integrator.p, integrator.t)
            du = Trixi.wrap_array(du_ode, integrator.p)
            u = Trixi.wrap_array(integrator.u, integrator.p)
            entropy_rate = Trixi.analyze(Trixi.entropy_timederivative, du, u, integrator.t, integrator.p)
            push!(data_entropy_rate, entropy_rate)
        end
    end

    filename = joinpath(DATA_DIR, "entropy_conservation_test_$(nameof(typeof(mesh)))$(ndims(mesh))D.h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["entropy_rate"] = data_entropy_rate
    end
    @info "Results saved" filename

    return filename
end

# Same as above but in 3D
function save_check_ec_3d(; initial_refinement_level = 5,
                            MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations3D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations3D)
        p = 1.0
        if iszero(x)
            v = zero(x)
        else
            v = -0.2 * normalize(x)
        end
        return prim2cons(SVector(v..., p), equations)
    end

    coordinates_min = (-2.0, -2.0, -2.0)
    coordinates_max = (+2.0, +2.0, +2.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level,
                        n_cells_max = 10^5,
                        periodicity = true)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level,
                        periodicity = true)
    end

    volume_flux = flux_thein_ranocha
    solver = DGSEM(polydeg = 3, surface_flux = volume_flux,
                   volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

    tspan = (0.0, 1.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 25)
    callbacks = CallbackSet(summary_callback, alive_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    data_t = Vector{Float64}()
    data_entropy_rate = Vector{Float64}()
    du_ode = similar(integrator.u)

    for _ in integrator
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            Trixi.rhs!(du_ode, integrator.u, integrator.p, integrator.t)
            du = Trixi.wrap_array(du_ode, integrator.p)
            u = Trixi.wrap_array(integrator.u, integrator.p)
            entropy_rate = Trixi.analyze(Trixi.entropy_timederivative, du, u, integrator.t, integrator.p)
            push!(data_entropy_rate, entropy_rate)
        end
    end

    filename = joinpath(DATA_DIR, "entropy_conservation_test_$(nameof(typeof(mesh)))$(ndims(mesh))D.h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["entropy_rate"] = data_entropy_rate
    end
    @info "Results saved" filename

    return filename
end


###########################################################
# Example 2 of Kunik, Kolb, Müller, and Thein (2024)
function save_self_similar_expansion_2d(; base_level = 3,
                                          med_level = 5,
                                          max_level = 9,
                                          amr_interval = 5,
                                          volume_flux = flux_thein_ranocha,
                                          MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations2D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations2D)
        p = 1.0
        if iszero(x)
            v = zero(x)
        else
            v = 0.5 * normalize(x)
        end
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else # P4estMesh, T8codeMesh
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc)
    end

    coordinates_min = (-2.0, -2.0)
    coordinates_max = (+2.0, +2.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^5,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 1.0,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 1.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 100)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    limiter! = PositivityPreservingLimiterZhangShu(thresholds = (1.0e-4,),
                                                   variables = (pressure,))

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    alg = SSPRK43(stage_limiter! = limiter!, step_limiter! = limiter!, thread = Trixi.True())
    integrator = init(ode, alg;
                      adaptive = true, controller = PIDController(0.55, -0.27, 0.05),
                      abstol = 1.0e-6, reltol = 1.0e-6, dt = 1.0e-8,
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 2, length = 1000) / sqrt(2)
    curve = vcat(x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    data_v = @. (v1 + v2) / sqrt(2)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            append!(data_v, @. (v1 + v2) / sqrt(2))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "self_similar_expansion_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end


# Same setup as Example 2 of Kunik, Kolb, Müller, and Thein (2024),
# but in 3D instead of 2D
function save_self_similar_expansion_3d(; base_level = 3,
                                          med_level = 5,
                                          max_level = 5,
                                          amr_interval = 5,
                                          volume_flux = flux_thein_ranocha,
                                          MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations3D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations3D)
        p = 1.0
        if iszero(x)
            v = zero(x)
        else
            v = 0.5 * normalize(x)
        end
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else # P4estMesh, T8codeMesh
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc,
                                   :z_neg => bc, :z_pos => bc)
    end

    coordinates_min = (-2.0, -2.0, -2.0)
    coordinates_max = (+2.0, +2.0, +2.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^6,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 1.0,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 1.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 25)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    limiter! = PositivityPreservingLimiterZhangShu(thresholds = (1.0e-8,),
                                                   variables = (pressure,))

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    alg = SSPRK43(stage_limiter! = limiter!, step_limiter! = limiter!, thread = Trixi.True())
    integrator = init(ode, alg;
                      adaptive = true, controller = PIDController(0.55, -0.27, 0.05),
                      abstol = 1.0e-8, reltol = 1.0e-8,
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 2, length = 1000) / sqrt(3)
    curve = vcat(x', x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    v3 = view(pd.data, :, 3)
    data_v = @. (v1 + v2 + v3) / sqrt(3)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            v3 = view(pd.data, :, 3)
            append!(data_v, @. (v1 + v2 + v3) / sqrt(3))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "self_similar_expansion_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end


###########################################################
# Example 3 of Kunik, Kolb, Müller, and Thein (2024)
function save_expansion_spherical_bubble_2d(; base_level = 3,
                                              med_level = 5,
                                              max_level = 9,
                                              amr_interval = 5,
                                              volume_flux = flux_thein_ranocha,
                                              MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations2D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations2D)
        if norm(x) <= 1
            p = 1.0
        else
            p = 0.1
        end
        v = zero(x)
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc)
    end

    coordinates_min = (-6.0, -6.0)
    coordinates_max = (+6.0, +6.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^5,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 0.5,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 6.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 100)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 6, length = 1000) / sqrt(2)
    curve = vcat(x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    data_v = @. (v1 + v2) / sqrt(2)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            append!(data_v, @. (v1 + v2) / sqrt(2))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "expansion_spherical_bubble_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end


# Same setup as Example 3 of Kunik, Kolb, Müller, and Thein (2024),
# but in 3D instead of 2D
function save_expansion_spherical_bubble_3d(; base_level = 3,
                                              med_level = 5,
                                              max_level = 7,
                                              amr_interval = 5,
                                              volume_flux = flux_thein_ranocha,
                                              MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations3D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations3D)
        if norm(x) <= 1
            p = 1.0
        else
            p = 0.1
        end
        v = zero(x)
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc,
                                   :z_neg => bc, :z_pos => bc)
    end

    coordinates_min = (-6.0, -6.0, -6.0)
    coordinates_max = (+6.0, +6.0, +6.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^6,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 0.5,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 6.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 25)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      abstol = 1.0e-6, reltol = 1.0e-6,
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 6, length = 1000) / sqrt(3)
    curve = vcat(x', x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    v3 = view(pd.data, :, 3)
    data_v = @. (v1 + v2 + v3) / sqrt(3)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            v3 = view(pd.data, :, 3)
            append!(data_v, @. (v1 + v2 + v3) / sqrt(3))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "expansion_spherical_bubble_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end


###########################################################
# Example 4 of Kunik, Kolb, Müller, and Thein (2024)
function save_collapse_spherical_bubble_2d(; base_level = 3,
                                             med_level = 5,
                                             max_level = 9,
                                             amr_interval = 5,
                                             volume_flux = flux_thein_ranocha,
                                             MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations2D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations2D)
        if norm(x) <= 1
            p = 0.1
        else
            p = 1.0
        end
        v = zero(x)
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc)
    end

    coordinates_min = (-6.0, -6.0)
    coordinates_max = (+6.0, +6.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^5,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 0.5,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 6.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 100)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 6, length = 1000) / sqrt(2)
    curve = vcat(x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    data_v = @. (v1 + v2) / sqrt(2)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            append!(data_v, @. (v1 + v2) / sqrt(2))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "collapse_spherical_bubble_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end


# Same setup as Example 4 of Kunik, Kolb, Müller, and Thein (2024),
# but in 3D instead of 2D
function save_collapse_spherical_bubble_3d(; base_level = 3,
                                             med_level = 5,
                                             max_level = 7,
                                             amr_interval = 5,
                                             volume_flux = flux_thein_ranocha,
                                             MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations3D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations3D)
        if norm(x) <= 1
            p = 0.1
        else
            p = 1.0
        end
        v = zero(x)
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc,
                                   :z_neg => bc, :z_pos => bc)
    end

    coordinates_min = (-6.0, -6.0, -6.0)
    coordinates_max = (+6.0, +6.0, +6.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^6,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 0.5,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 6.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 25)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 6, length = 1000) / sqrt(3)
    curve = vcat(x', x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    v3 = view(pd.data, :, 3)
    data_v = @. (v1 + v2 + v3) / sqrt(3)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            v3 = view(pd.data, :, 3)
            append!(data_v, @. (v1 + v2 + v3) / sqrt(3))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "collapse_spherical_bubble_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end



###########################################################
# Example 5 of Kunik, Kolb, Müller, and Thein (2024)
function save_sine_velocity_2d(; base_level = 3,
                                 med_level = 5,
                                 max_level = 9,
                                 amr_interval = 5,
                                 volume_flux = flux_thein_ranocha,
                                 MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations2D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations2D)
        p = 1.0
        r = norm(x)
        if r < 1
            v = 2 * π * sinc(2 * r) * x
        else
            v = zero(x)
        end
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc)
    end

    coordinates_min = (-5.0, -5.0)
    coordinates_max = (+5.0, +5.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^5,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 0.5,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 6.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 100)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 5, length = 1000) / sqrt(2)
    curve = vcat(x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    data_v = @. (v1 + v2) / sqrt(2)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            append!(data_v, @. (v1 + v2) / sqrt(2))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "sine_velocity_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end


# Same setup as Example 5 of Kunik, Kolb, Müller, and Thein (2024),
# but in 3D instead of 2D
function save_sine_velocity_3d(; base_level = 3,
                                 med_level = 5,
                                 max_level = 7,
                                 amr_interval = 5,
                                 volume_flux = flux_thein_ranocha,
                                 MeshType = P4estMesh)
    equations = UltraRelativisticEulerEquations3D()

    function initial_condition(x, t, equations::UltraRelativisticEulerEquations3D)
        p = 1.0
        r = norm(x)
        if r < 1
            v = 2 * π * sinc(2 * r) * x
        else
            v = zero(x)
        end
        return prim2cons(SVector(v..., p), equations)
    end

    if MeshType === TreeMesh
        boundary_conditions = BoundaryConditionDirichlet(initial_condition)
    else
        bc = BoundaryConditionDirichlet(initial_condition)
        boundary_conditions = Dict(:x_neg => bc, :x_pos => bc,
                                   :y_neg => bc, :y_pos => bc,
                                   :z_neg => bc, :z_pos => bc)
    end

    coordinates_min = (-5.0, -5.0, -5.0)
    coordinates_max = (+5.0, +5.0, +5.0)
    if MeshType === TreeMesh
        mesh = TreeMesh(coordinates_min, coordinates_max;
                        initial_refinement_level = base_level,
                        n_cells_max = 10^6,
                        periodicity = false)
    else # P4estMesh, T8codeMesh
        mesh = MeshType((1, 1, 1); polydeg = 1,
                        coordinates_min, coordinates_max,
                        initial_refinement_level = base_level,
                        periodicity = false)
    end

    surface_flux = flux_lax_friedrichs
    basis = LobattoLegendreBasis(3)
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 0.5,
                                             alpha_min = 0.001,
                                             alpha_smooth = true,
                                             variable = pressure_velocity)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
    solver = DGSEM(basis, surface_flux, volume_integral)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                        boundary_conditions)

    tspan = (0.0, 6.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()
    alive_callback = AliveCallback(alive_interval = 25)

    amr_indicator = IndicatorLöhner(semi, variable = pressure_velocity)
    amr_controller = ControllerThreeLevel(semi, amr_indicator;
                                          base_level,
                                          med_level, med_threshold = 0.05,
                                          max_level, max_threshold = 0.1)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = amr_interval,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    callbacks = CallbackSet(summary_callback, alive_callback, amr_callback)

    # SSPRK43 with optimized controller of Ranocha, Dalcin, Parsani,
    # and Ketcheson (2021)
    integrator = init(ode, SSPRK43(thread = Trixi.True());
                      controller = PIDController(0.55, -0.27, 0.05),
                      abstol = 1.0e-6, reltol = 1.0e-6,
                      callback = callbacks, ode_default_options()...)

    # Prepare plotting
    x = range(0, 5, length = 1000) / sqrt(3)
    curve = vcat(x', x', x')
    pd = PlotData1D(integrator.u, semi; curve)
    data_r = pd.x
    data_t = [integrator.t]
    v1 = view(pd.data, :, 1)
    v2 = view(pd.data, :, 2)
    v3 = view(pd.data, :, 3)
    data_v = @. (v1 + v2 + v3) / sqrt(3)
    data_p = pd.data[:, end]

    tstops = range(ode.tspan..., length = length(x))
    for _ in TimeChoiceIterator(integrator, tstops)
        push!(data_t, integrator.t)
        Trixi.@trixi_timeit Trixi.timer() "Prepare data for plotting" begin
            pd = PlotData1D(integrator.u, semi; curve)
            @assert data_r ≈ pd.x
            v1 = view(pd.data, :, 1)
            v2 = view(pd.data, :, 2)
            v3 = view(pd.data, :, 3)
            append!(data_v, @. (v1 + v2 + v3) / sqrt(3))
            append!(data_p, pd.data[:, end])
        end
    end

    summary_callback()

    filename = joinpath(DATA_DIR, "sine_velocity_$(nameof(typeof(mesh)))$(ndims(mesh))D_$(base_level)_$(med_level)_$(max_level).h5")
    h5open(filename, "w") do file
        file["time"] = data_t
        file["radius"] = data_r
        file["velocity"] = reshape(data_v, length(data_r), length(data_t))
        file["pressure"] = reshape(data_p, length(data_r), length(data_t))
    end
    @info "Results saved" filename

    return filename
end
