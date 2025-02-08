# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

module ITPRootFinder

const GOLDEN_RATIO = Base.MathConstants.golden

include("main.jl")

public Config, find_root!

end # module ITPRootFinder
