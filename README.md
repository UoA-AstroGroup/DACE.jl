# DACE.jl

DACE.jl provides differential algebra in Julia using direct bindings to the
[upstream DACE C engine](https://github.com/dacelib/dace).
Polynomial arithmetic runs in DACE; evaluation, map inversion, arrays, eigenpairs
and moments are implemented in Julia. CxxWrap and Eigen are no longer required.

This branch builds the native engine from a pinned upstream commit. From this
checkout, with Julia 1.10 or later, Git, CMake and a C/C++ compiler installed:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.build("DACE")
Pkg.test()
```

On Windows, use MinGW-w64, or the supported fallback: WSL Ubuntu with CMake and
`gcc-mingw-w64-x86-64-posix` / `g++-mingw-w64-x86-64-posix` installed.
A compiler-free installation will require a new C-only binary artifact; this
branch does not use the existing CxxWrap-based DACE_jll.

```julia
using DACE
DACE.init(10, 1)
x = DA(1, 1.0)                  # Independent variable; DA(1) is constant one.
y = sin(x)
DACE.evaluate(y, [0.1])
DACE.getCoefficient(y, [3])     # -1/6
```

Polynomial evaluation is now named `DACE.evaluate` (previously `DACE.eval`).
`evalScalar` is retained. See the [binding guide](docs/src/tutorials/direct-bindings.md)
for API coverage, memory ownership, concurrency and performance, and
[examples](examples) for polynomial inversion and ODE integration.
