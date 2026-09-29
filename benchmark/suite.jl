# Separate processes are required: julia --project=... suite.jl direct|cxxwrap [results.toml]
# Set DACE_IOD_EXAMPLES to the folder containing tutorials 08 and 09.
using Printf, TOML, LinearAlgebra, SHA
const backend = ARGS[1]
backend in ("direct", "cxxwrap") || error("Choose direct or cxxwrap")
import DACE
const D = DACE
const DA = D.DA
LinearAlgebra.BLAS.set_num_threads(1)

wrap(v) = backend == "direct" ? v : D.AlgebraicVector(v)
evaluate(p, x) = backend == "direct" ? D.evaluate(p, x) : D.eval(p, x)
if backend == "cxxwrap"
    # Select the native overload explicitly: the current package has an ambiguity
    # between compiledDA/AbstractVector and the wrapped AlgebraicVector overload.
    evaluate(p::D.compiledDA, x) = D.eval(D.CxxWrap.ConstCxxRef(p), x)
end
coeffnorm(p) = D.norm(p, 0)
terms(c) = backend == "direct" ? length(c.levels) : D.getTerms(c)

function repeat_job(f::F, n) where F
    checksum = 0.0
    for _ in 1:n
        checksum += f()
    end
    isfinite(checksum) || error("Nonfinite benchmark checksum")
    checksum
end

const rows = Dict{String,Any}[]
function measure(label, f::F; limit=10_000) where F
    repeat_job(f, 2) # Compile, then warm.
    GC.gc(true)
    pilot = @elapsed repeat_job(f, 5)
    n = clamp(ceil(Int, 0.03 / (pilot/5)), 1, limit)
    times, bytes, cleanup = Float64[], Float64[], Float64[]
    for _ in 1:9
        GC.gc(true)
        stats = @timed repeat_job(f, n)
        push!(times, stats.time/n)
        push!(bytes, stats.bytes/n)
        push!(cleanup, (@elapsed GC.gc(true))/n)
    end
    med = sort(times)[5]
    println(rpad(label, 30), @sprintf(" %10.3f us/op", 1e6*med))
    push!(rows, Dict("name"=>label, "seconds"=>times, "median_seconds"=>med,
        "julia_bytes_per_op"=>sort(bytes)[5], "cleanup_seconds"=>cleanup, "iterations"=>n))
end

function numeric_job(c, nv, nout)
    if backend == "direct"
        a, out, work = fill(1e-3, nv), zeros(nout), zeros(c.depth+1)
        return () -> (D.evaluate!(out, c, a, work); out[1])
    else
        a, out = wrap(fill(1e-3, nv)), wrap(zeros(nout))
        return () -> (D.eval(c, a, out); out[1])
    end
end

function algebra(no, nv)
    D.init(no, nv)
    D.setEps(1e-24)
    x = [DA(i, 1.0) for i in 1:nv]
    p, q = 0.1 + sum(x)/nv, 0.2 + sum(i*x[i] for i in 1:nv)/nv
    kernel() = sin(p*q) + exp(p/10)/(2 + q*q)
    @assert D.cons(kernel()) ≈ sin(0.02) + exp(0.01)/2.04
    prefix = "o$(no)v$(nv) "
    measure(prefix*"arithmetic", () -> D.cons(kernel()))
    values = wrap([sin(p*q), exp(p/10)/(2 + q*q)])
    c = D.compile(values)
    measure(prefix*"compile map", () -> Float64(terms(D.compile(values))); limit=2_000)
    evaljob = numeric_job(c, nv, 2)
    @assert evaljob() ≈ sin((0.1+1e-3)*(0.2+sum(1:nv)*1e-3/nv)) atol=1e-13
    measure(prefix*"numeric evaluation", evaljob; limit=500_000)
    args = wrap([x[i] + 0.01x[mod1(i+1,nv)]^2 for i in 1:nv])
    measure(prefix*"DA composition", () -> D.cons(evaluate(c, args)[1]); limit=500)
    f = wrap([x[i] + 0.02*(sum(x)/nv)^2 + 0.01*x[i]*x[mod1(i+1,nv)] for i in 1:nv])
    inverse = D.invert(f)
    residual = collect(evaluate(f, inverse)) - x
    @assert maximum(coeffnorm, residual) < 1e-12
    measure(prefix*"map inversion", () -> D.cons(D.invert(f)[1]); limit=200)
end

function load_tutorial(name, filename)
    path = joinpath(ENV["DACE_IOD_EXAMPLES"], filename)
    source = read(path, String)
    source = replace(source, r"\nmain\(\)\s*$" => "\n")
    if backend == "direct"
        source = replace(source, "DACE.eval(" => "DACE.evaluate(")
    end
    mod = Module(name)
    Base.include_string(mod, source, path)
    return mod
end

if haskey(ENV, "DACE_IOD_EXAMPLES")
    const Kepler = load_tutorial(:Kepler, "08_dace_implicit_kepler.jl")
    const IOD = load_tutorial(:IOD, "09_daiod_dace.jl")
end

function workloads()
    # Validate each complete tutorial first. Suppress console I/O for timings.
    redirect_stdout(devnull) do
        Kepler.main()
        IOD.main()
    end
    times, mu = [-120.0, 0.0, 150.0], 398600.4418
    angles, sites = IOD.observations(times)
    seed = only(IOD.gauss(angles, sites, times, mu)).rho
    nominal = IOD.match_ranges(seed, angles, sites, times, mu)
    # Keep numerical assertions in the workloads; exclude scalar setup from map timing.
    redirect_stdout(devnull) do
        measure("Kepler implicit map", () -> D.cons(Kepler.implicit_kepler(1.3, 0.3, 0.01)); limit=100)
        measure("DAIOD map construction", () -> begin
            orbit = IOD.daiod(nominal.rho, angles, sites, times, mu)
            D.cons(orbit.state[1])
        end; limit=10)
    end
    # measure printed into devnull above; show its stored results here.
    for row in rows[end-1:end]
        println(rpad(row["name"], 30), @sprintf(" %10.3f us/op", 1e6*row["median_seconds"]))
    end
end

println("Backend: ", backend, "; Julia ", VERSION, "; ", Sys.MACHINE,
    "; threads=", Threads.nthreads(), "; BLAS threads=", BLAS.get_num_threads())
algebra(6, 2)
algebra(3, 9)
haskey(ENV, "DACE_IOD_EXAMPLES") && workloads()
if length(ARGS) == 2
    native = backend == "direct" ? D.lib : D.DACE_jll.libdace
    open(ARGS[2], "w") do io
        TOML.print(io, Dict("backend"=>backend, "julia"=>string(VERSION),
            "machine"=>Sys.MACHINE, "rows"=>rows, "native_library"=>native,
            "native_sha256"=>bytes2hex(open(sha256, native)),
            "julia_threads"=>Threads.nthreads(), "blas_threads"=>BLAS.get_num_threads()))
    end
end
