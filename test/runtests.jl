# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

using Test, Random
import ITPRootFinder as ITP

@testset "cbrt(x) any finite initial start" begin

    # root at x == 0.
    function eval_func2(x)
        return cbrt(x)
    end

    struct FuncCallable{FT}
        f::FT
    end
    function (A::FuncCallable)(_, x, _)
        return A.f(x)
    end

    T = Float64
    abs_tol = T(1e-3)
    Z = 10_000
    rng = Random.Xoshiro(0)

    config = ITP.Config(T; f_tol=abs_tol, x_tol=abs_tol, max_iters=typemax(Int))

    N_tests = 1000
    for _ = 1:N_tests
        lb = -rand(rng, T)*2Z
        ub = rand(rng, T)*Z
        r, status = ITP.find_root!(nothing, FuncCallable(eval_func2), nothing, lb, ub, config)

        @test isapprox(r,0,atol=abs_tol)
    end

end
