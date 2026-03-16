# DACE.jl

This is the documentation page for [DACE.jl](https://github.com/UoA-AstroGroup/DACE.jl),
a Julia interface to the [DACE library](https://github.com/dacelib/dace).

## Getting started

DACE.jl can be installed using the Julia package manager. From the Julia REPL,
type `]` to enter the Pkg REPL mode and run

```
pkg> add DACE
```

## Notes about the interface

- The Julia interface is built using [CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl)
- The C++ source for the interface is currently in [this fork](https://github.com/UoA-AstroGroup/dace/tree/julia-interface/interfaces/julia)
  of the DACE library
- The C++ code gets built by [BinaryBuilder](https://docs.binarybuilder.org/stable/) and released into the Julia registry
  - the build recipe is located [here](https://github.com/JuliaPackaging/Yggdrasil/tree/master/D/DACE)
- The Julia component of the interface is currently [here](https://github.com/UoA-AstroGroup/DACE.jl)
  - the DACE library has been released to the Julia General registry

The above may change, in particular we hope to:

- merge the forked DACE library back to [upstream](https://github.com/dacelib/dace), if possible
- move the DACE.jl package under the [dacelib](https://github.com/dacelib) organisation on GitHub, if possible

## Differences compared to the C++ interface

### DA constructor

The `DA` constructor is different in the Julia interface compared to C++. In C++ if you create a
`DA` object with a single argument the behaviour is different depending on whether the argument is
an integer or double. For example, in C++:

- `DA(1)` creates a DA object representing 1.0 times the independent variable number 1
- `DA(1.0)` creates a DA object with the constant part equal to 1.0
- `DA(1, 1.0)` or `DA(1, 1)` both create a DA representing 1.0 times the independent variable number 1 (same as first bullet point above)

In Julia with CxxWrap, there is an issue detecting the difference between integer and double when passing a single value on Windows,
resulting in huge numbers of warnings being printed to screen. Therefore, in the Julia interface we have:

- `DA(1)` and `DA(1.0)` both create a DA object with the constant part equal to 1.0
- `DA(1, 1.0)` or `DA(1, 1)` both create a DA representing 1.0 times the independent variable number 1
