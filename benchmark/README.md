# Performance comparison

`julia --project=. benchmark/suite.jl direct results/direct.toml` benchmarks scalar
arithmetic, map compilation, numeric evaluation, DA composition and inversion.
Each kernel is warmed, calibrated and timed in nine batches. Natural GC during
each batch is included. Full GC before and after a batch is recorded separately.
Numeric evaluation reuses its output and workspace.

To compare with CxxWrap, use a separate Julia environment containing the old
DACE.jl implementation and its matching binary. Run the same `suite.jl` there
with `cxxwrap` instead of `direct`. The harness explicitly selects the native
CxxWrap evaluation overload where the old package has ambiguous methods.

For paired Linux/WSL measurements, set `DACE_CXX_PROJECT_110` and
`DACE_CXX_PROJECT_113` to those environments, then run:

```sh
julia benchmark/run_linux.jl
julia benchmark/summarize.jl
```

The runner uses Juliaup, CPU 2, one Julia thread and one BLAS thread. It alternates
backends and reverses their order in the second round. `DACE_JULIA` and
`DACE_BENCH_CPU` override the launcher and CPU. `DACE_DIRECT_PROJECT_110` and
`DACE_DIRECT_PROJECT_113` can select separate resolved environments for this
checkout on each Julia version. Run on an otherwise idle machine.
Raw samples and native library hashes are stored under the ignored `results/paired`.

Optionally set `DACE_IOD_EXAMPLES` to DSTOrbitDetermination.jl's `examples/iod`
folder to include the implicit Kepler and DAIOD tutorials. The harness executes
their numerical validations before timing. DAIOD timing measures construction of
the angle-to-state map, excluding the initial scalar Gauss/range solution.

Julia allocation counts are not directly comparable: this implementation exposes
its coefficient buffers to Julia's collector; the CxxWrap backend allocates those
buffers in native memory. Timings compare the complete implementations, including
their native compiler and threading configuration, rather than isolating FFI cost.
