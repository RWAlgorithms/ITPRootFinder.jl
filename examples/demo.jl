# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

using Random, LinearAlgebra, Statistics
import PythonPlot as PLT

using Revise
import ITPRootFinder as ITP

PLT.close("all")
fig_num = 1

# root at x == +/-1.
function eval_func1(x)
    return x^20 − 1
end

# root at x == 0.
function eval_func2(x)
    return cbrt(x)
end

# root around x == −1.76929
function eval_func3(x)
    return x^3 - 2 * x + 2
end

struct FuncCallable{FT}
    f::FT
end
function (A::FuncCallable)(_, x, _)
    return A.f(x)
end

T = Float64
abs_tol = T(1e-3)

# # Demo
config = ITP.Config(T; f_tol=abs_tol, x_tol=abs_tol, max_iters=typemax(Int))

lb = T(0.1)
ub = T(5)
r, status = ITP.find_root!(nothing, FuncCallable(eval_func1), nothing, lb, ub, config)
@show lb, ub, r, status

# # Visualize.
ts = LinRange(lb, ub, 10_000)
f = tt -> eval_func1(tt)
f_ts = f.(ts)

PLT.figure(fig_num)
fig_num += 1
PLT.plot(ts, f_ts)
PLT.title("Function 1")

# # Second test func.

lb = T(-10)
ub = T(2)
r, status = ITP.find_root!(nothing, FuncCallable(eval_func2), nothing, lb, ub, config)
@show lb, ub, r, status
println("The above status was false because x_tol was reached before f_tol.")

# # Visualize.
ts = LinRange(lb, ub, 10_000)
f = tt -> eval_func2(tt)
f_ts = f.(ts)

PLT.figure(fig_num)
fig_num += 1
PLT.plot(ts, f_ts)
PLT.title("Function 2")

# # Third test func.

lb = T(-100)
ub = T(5)
r, status = ITP.find_root!(nothing, FuncCallable(eval_func3), nothing, lb, ub, config)
@show lb, ub, r, status

# # Visualize.
ts = LinRange(lb, ub, 10_000)
f = tt -> eval_func3(tt)
f_ts = f.(ts)

PLT.figure(fig_num)
fig_num += 1
PLT.plot(ts, f_ts)
PLT.title("Function 3")


using BenchmarkTools
@btime ITP.find_root!(nothing, FuncCallable($eval_func3), nothing, $lb, $ub, $config)

nothing
