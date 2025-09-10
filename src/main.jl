# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

"""
    Config(
        ::Type{T};
        f_tol::T=T(1e-4),
        x_tol::T=T(1e-6),
        k1::T=T(0.1),
        k2::T=T(0.98 * (1 + GOLDEN_RATIO)), # see equation 24.
        n0::Int=0,
        max_iters::Int=typemax(Int),
    ) where {T<:AbstractFloat}
Usually the ITP algorithm terminates either `x_tol` or `f_tol` is reached, but `max_iters` is an additional early stopping condition. Please consult the work of *Oliveira and Takahashi, 2020, DOI: 10.1145/3423597* for the interpretation of `k1`, `k2`, `n0`.
"""
struct Config{T <: AbstractFloat}
    f_tol::T
    x_tol::T
    k1::T
    k2::T
    n0::Int
    max_iters::Int # guard.

    function Config(
            ::Type{T};
            f_tol::T = T(1.0e-4),
            x_tol::T = T(1.0e-6),
            k1::T = T(0.1),
            k2::T = T(0.98 * (1 + GOLDEN_RATIO)), # see equation 24.
            n0::Int = 0,
            max_iters::Int = typemax(Int),
        ) where {T <: AbstractFloat}

        k1 > zero(T) || error("The following must be true: k1 > 0")
        one(T) <= k2 < one(T) + GOLDEN_RATIO || error("The following must be true: 1 <= k2 < 1 + golden ratio")
        n0 >= 0 || error("The following must be true: n0 >= 0")
        x_tol > 0 || error("The following must be true: x_tol > 0")
        f_tol > 0 || error("The following must be true: f_tol > 0")

        return new{T}(f_tol, x_tol, k1, k2, n0, max_iters)
    end
end


"""
    function find_root!(
        f_buffer, # potentially mutates
        f!, # a callable
        f_params,
        lb::T,
        ub::T,
        config::Config{T},
    ) where {T<:AbstractFloat}

Finds **one** root in `[lb,ub]`, to given `abs(f(x_next)-f(x_current)) < f_tol` or `abs(x_next-x_current)<=x_tol` stopping conditions. The interval `[lb,ub]` must be such that `f(lb) * f(ub) < 0` and `lb < ub`.

- Based on *Oliveira and Takahashi, 2020, DOI: 10.1145/3423597*.
- Assumes `f(lb) < 0 < f(ub)`.
- f! should be a callable of the form `f!(f_buffer, x, f_params)` that returns a scalar of type `T`. `x` is the scalar input, `f_buffer` is potentially mutated when computing `f!`, and `f_params` remains unchanged as one calls `f` with different values for `x`.

Returns `(root_sol, status)`, where `root_sol` is of type `T`, and `status` is a `Bool`.
- `status==true` if `f_tol` was reached.
- If `x_tol` was reached before `f_tol`, then `status==false` and `isfinite(root_sol)==true`.
- If `status==false` and `isfinite(root_sol)==false`, then the assumptions `f_lb * f_ub < 0` and  `lb < ub` are violated. This means `find_root!` cannot find a root.
"""
function find_root!(
        f_buffer, # potentially mutates
        f!, # a callable
        f_params,
        lb::T,
        ub::T,
        config::Config{T},
    ) where {T <: AbstractFloat}

    lb < ub || error("lb must be smaller than ub.")

    # check if we have a valid problem to solve. See equation 2.
    f_lb = f!(f_buffer, lb, f_params)
    f_ub = f!(f_buffer, ub, f_params)
    if !(f_lb * f_ub < 0) || lb > ub
        # cannot do root finding on this interval.
        return T(NaN), false
    end

    # set up.
    f_tol, x_tol, n0, k1, k2, max_iters = config.f_tol, config.x_tol, config.n0, config.k1, config.k2, config.max_iters

    ϵ = x_tol / 2
    n_bin = ceil(Int, log2((ub - lb) / x_tol)) # Theorem 1.1.
    n_max = n_bin + n0

    # # bracketing algorithm.
    k = 0
    a = lb
    b = ub
    f_a = f_lb
    f_b = f_ub
    C = 2^(n_max + 1)

    while abs(b - a) > x_tol && k < max_iters

        # interpolation.
        x_RF = (a * f_b - b * f_a) / (f_b - f_a) # equation 5.

        # truncation, x_t.
        x_bin = (a + b) / 2 # equation 4

        σ = sign(x_bin - x_RF)
        δ = k1 * abs(b - a)^k2

        x_t = x_bin
        if δ <= abs(x_bin - x_RF)
            x_t = x_RF + σ * δ
        end

        # minimax radius and interval.
        #r_k = ϵ*2^(n_bin-k) - (b-a)/2
        C = C / 2
        r_k = ϵ * C - (b - a) / 2

        #interval_k = [x_bin - r_k, x_bin + r_k]

        # projection onto interval_k.
        x_ITP = x_t
        if abs(x_t - x_bin) > r_k
            x_ITP = x_bin - σ * r_k
        end

        # choose a new position in the interval, and see if it is numerically close to being a root.
        # w = x_ITP # unguarded.
        w = clamp(x_ITP, lb, ub) # guarded.
        #w = (a + b)/2 # uncomment to overide ITP with bisection.

        f_w = f!(f_buffer, w, f_params) #f(w)

        if abs(f_w) < f_tol
            return w, true
        end

        # Set up for the next iteration: bracket steps.
        if f_w * f_a > 0
            a = w
            f_a = f_w
            # b = b
        elseif f_w * f_b > 0
            #a = a
            b = w
            f_b = f_w
        else
            a = w
            b = w
        end

        # @show k # timing diagnostics
        k += 1 # first iteration is defined as k == 0. See the paragraph after equation 16.
    end

    return clamp((a + b) / 2, lb, ub), false
end
