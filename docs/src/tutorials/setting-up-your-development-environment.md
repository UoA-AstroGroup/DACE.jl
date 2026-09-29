# Setting up your development environment

Use Julia 1.10 or later. Install Git, CMake, and C/C++ compilers (for example,
`build-essential cmake git` on Ubuntu, or Xcode command-line tools and CMake on
macOS). The upstream CMake project detects both compilers, but only its C core
is linked into this package.

From a local checkout of DACE.jl:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.build("DACE")
Pkg.test(; julia_args=["--threads=4"])
```

On Windows, put CMake and MinGW-w64 on PATH. If they are unavailable, the build
uses WSL Ubuntu; install `cmake`, `gcc-mingw-w64-x86-64-posix` and
`g++-mingw-w64-x86-64-posix` in that distribution. Run the build from Windows Julia
to produce a Windows DLL, or from Linux Julia to produce a Linux shared library.

`deps/build.jl` checks out a pinned upstream revision in `deps/source-upstream`,
builds with CMake and installs under `deps/usr`. An existing checkout at the same
revision can be passed to `julia deps/build.jl /path/to/dace`. Build products and
cloned sources are ignored by Git. No Julia interface branch of DACE is needed.

To build documentation locally, develop this checkout in the docs environment:

```julia
Pkg.activate("docs")
Pkg.develop(path=".")
Pkg.instantiate()
include("docs/make.jl")
```

The documentation build executes the sine, inversion and ODE examples. Deployment
is only attempted in CI. Before a public release, package the C-only build as a
binary artifact so users do not need a compiler.
