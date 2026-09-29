# DACE.jl

[DACE.jl](https://github.com/UoA-AstroGroup/DACE.jl) binds directly to the
[upstream DACE C engine](https://github.com/dacelib/dace). Julia implements the
higher-level array, evaluation, inversion, eigenpair and moment operations.

## Getting started

This branch requires Julia 1.10 or later and builds a pinned native library from
source. Follow [Setting up your development environment](tutorials/setting-up-your-development-environment.md).

```julia
using DACE
DACE.init(6, 2)
x, y = DACE.identity()
p = sin(x) * exp(y)
DACE.evaluate(p, [0.1, 0.2])
DACE.getCoefficient(p, [1, 1])
```

`DA(c)` creates a constant, including when `c` is an integer. `DA(i, c)` creates
`c` times independent variable `i`; index zero creates a constant.

Polynomial evaluation uses `DACE.evaluate`; `evalScalar` remains available.
See [Direct bindings](tutorials/direct-bindings.md) for compatibility and ownership details.
