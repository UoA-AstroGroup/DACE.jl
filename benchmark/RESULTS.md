# Local validation and benchmarks

Validated on 2026-09-29 with Julia 1.10.12 and 1.13.1:

- All **950 tests pass** on Windows x86-64 and Linux x86-64 (WSL Ubuntu), with four Julia threads. This includes the original mathematical tests and 169 new binding checks.
- Native builds succeed for both operating systems using unmodified upstream DACE. Each binary works with both Julia versions.
- The documentation builds locally, executing sine, polynomial inversion and the orbital ODE example. The ODE's nominal endpoint differs from the ordinary numeric solution by `5.12e-15`. The DiffEqBase extension also passes its value-extraction check. Local WSL documentation validation disables Git source links because the worktree is owned by Windows.
- The Kepler and DAIOD tutorial validations pass for both backends on both Julia versions before benchmarking.
- Numeric `evaluate!` allocates zero Julia bytes with reusable buffers and all coordinates supplied.

Julia-managed coefficient buffers replace native finalizers. Reserving less space for constants and sparse arithmetic reduces DAIOD's measured Julia allocations from **16.70 MB to 7.12 MB per map** on 1.13 (57%). This comparison is within the new implementation, before and after that optimization. CxxWrap's native coefficient allocations are invisible to Julia's byte counter and must not be compared directly with these figures.

## Timing results

The tables report medians of 18 batch measurements: nine batches in each of two fresh processes, with backend order reversed for the second round. Runs use an Intel i7-12700, CPU affinity 2, one Julia thread and one BLAS thread. Timings include natural GC within batches, and exclude compilation, explicit full GC between batches and console output. `o6v2` denotes order 6 in two variables; `o3v9` denotes order 3 in nine variables. DAIOD timing excludes the initial scalar range solution.

**These runs were noisy.** DAIOD process medians were 3.40/2.93 ms (direct) versus 3.47/3.93 ms (CxxWrap) on 1.13, and 5.23/3.93 ms versus 7.16/5.33 ms on 1.10. Some synthetic cases varied by approximately 2x. Do not interpret small differences as established speedups. The pooled results favor direct bindings for the larger composition/inversion workloads and DAIOD, but small-map inversion and compilation regress on 1.10, and some numeric evaluation cases are slower.

The CxxWrap baseline is DACE.jl commit `d390ea198b2bf1aa127f599a57376cf118d2c361`, with the local DACE_jll 0.7.2 build, GCC 12.1 and pthread support. Direct bindings use upstream DACE `368f6020342b0223ef1c05d50602bc54b7a8fe1f`, GCC 11.4, pthread support disabled and a Julia lock around engine calls. This compares implementations and build configurations, rather than isolating FFI overhead. No Windows timing claims are made.

Reproduction instructions are in [README.md](README.md). Raw samples and native binary hashes are retained in the local ignored `results/paired` directory.

Julia 1.13.1
| Workload | CxxWrap (us) | Direct (us) | Direct/CxxWrap |
| --- | ---: | ---: | ---: |
| o6v2 arithmetic | 4.243 | 2.675 | 0.631 |
| o6v2 compile map | 1.631 | 1.207 | 0.740 |
| o6v2 numeric evaluation | 0.136 | 0.108 | 0.790 |
| o6v2 DA composition | 6.577 | 3.411 | 0.519 |
| o6v2 map inversion | 17.359 | 11.900 | 0.686 |
| o3v9 arithmetic | 12.076 | 10.317 | 0.854 |
| o3v9 compile map | 16.562 | 16.119 | 0.973 |
| o3v9 numeric evaluation | 0.564 | 0.817 | 1.450 |
| o3v9 DA composition | 129.183 | 113.832 | 0.881 |
| o3v9 map inversion | 780.719 | 542.630 | 0.695 |
| Kepler implicit map | 67.015 | 46.619 | 0.696 |
| DAIOD map construction | 3763.495 | 3249.430 | 0.863 |

Julia 1.10.12
| Workload | CxxWrap (us) | Direct (us) | Direct/CxxWrap |
| --- | ---: | ---: | ---: |
| o6v2 arithmetic | 5.027 | 4.179 | 0.831 |
| o6v2 compile map | 1.524 | 2.251 | 1.477 |
| o6v2 numeric evaluation | 0.116 | 0.136 | 1.169 |
| o6v2 DA composition | 7.946 | 6.456 | 0.812 |
| o6v2 map inversion | 19.227 | 26.959 | 1.402 |
| o3v9 arithmetic | 14.818 | 12.778 | 0.862 |
| o3v9 compile map | 19.775 | 21.531 | 1.089 |
| o3v9 numeric evaluation | 0.618 | 0.632 | 1.023 |
| o3v9 DA composition | 143.611 | 98.084 | 0.683 |
| o3v9 map inversion | 997.700 | 516.154 | 0.517 |
| Kepler implicit map | 93.625 | 90.312 | 0.965 |
| DAIOD map construction | 5578.630 | 4489.762 | 0.805 |
