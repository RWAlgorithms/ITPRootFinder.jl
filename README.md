# ITPRootFinder.jl
Consider a given callable `f!(f_buffer, x, params)` and an interval `[lb,ub]` such that `f(lb) * f(ub) < 0` and `lb < ub`. This implies `lb` cannot be a root, and `ub` cannot be a root.

`f` should be a callable of the form `f!(f_buffer, x, f_params)` that returns a scalar of type `T<:AbstractFloat`. `x` is the scalar input, `f_buffer` is potentially mutated when computing `f!`, and `f_params` remains unchanged as one calls `f` with different values for `x`. This separation of mutating (optional) input `f_buffer`, the actual function input `x`, and constants (with respect to the lifetime of `f`) `f_params`, is to allow users to pass in callables that attempt to avoid issues from the *world-age performance-related problem* in Julia. For a pure function, `f_buffer` and `f_params` can be set to anything, e.g., `nothing`.

The `find_root!` function in this package finds one root in `[lb,ub]`, to given `abs(f(x_next)-f(x_current)) < f_tol` or `abs(x_next-x_current)<=x_tol` stopping conditions.

For example:
```{julia}
import ITPRootFinder as ITP

# roots are at x == 1, x==-2.
function eval_func1(x)
    return x^20 − 1
end

struct FuncCallable{FT}
    f::FT
end
function (A::FuncCallable)(_, x, _)
    return A.f(x)
end

T = Float64
abs_tol = T(1e-3)

# bounds.
lb = T(0.1)
ub = T(5)

# configuration.
config = ITP.Config(T; f_tol=abs_tol, x_tol=abs_tol, max_iters=typemax(Int))

r, status = ITP.find_root!(nothing, FuncCallable(eval_func1), nothing, lb, ub, config)
@show lb, ub, r, status
```
The REPL should display something like this:
```
(lb, ub, r, status) = (0.1, 5.0, 0.9999902260018678, true)
```

# Citation
If you find this software useful, please cite this repository and code via the *Cite This Repository* button on the GitHub webpage of this repository.

# License
This project is licensed under the GPL V3.0 license; see the `LICENSE` file for details. Individual source files may contain the following tag instead of the full license text:
```
SPDX-License-Identifier: GPL-3.0-only
```

Using SPDX enables machine processing of license information based on the SPDX License Identifiers and makes it easier for developers to see at a glance which license they are dealing with.
