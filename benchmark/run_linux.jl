# Runs sequentially, alternating backends, pinned to one CPU.
# Set DACE_IOD_EXAMPLES and DACE_CXX_PROJECT_110 / DACE_CXX_PROJECT_113.
Sys.islinux() || error("Run in Linux/WSL")
root = normpath(joinpath(@__DIR__, ".."))
output = joinpath(root, "results", "paired")
mkpath(output)
juliaup = get(ENV, "DACE_JULIA", joinpath(homedir(), ".juliaup/bin/julia"))
cpu = get(ENV, "DACE_BENCH_CPU", "2")
for version in ("1.13", "1.10"), round in 1:2
    # Reverse the order in the second round to reduce systematic ordering bias.
    backends = round == 1 ? ("cxxwrap", "direct") : ("direct", "cxxwrap")
    for backend in backends
        tag = replace(version, "." => "")
        project = backend == "direct" ? get(ENV, "DACE_DIRECT_PROJECT_"*tag, root) : ENV["DACE_CXX_PROJECT_"*tag]
        stem = joinpath(output, "$(backend)-$(tag)-$(round)")
        println("Running Julia ", version, " / ", backend, " / round ", round)
        flush(stdout)
        cmd = `taskset -c $cpu $juliaup +$version --startup-file=no --threads=1 --project=$project
            $(@__DIR__)/suite.jl $backend $stem.toml`
        open(stem*".log", "w") do io
            run(pipeline(cmd; stdout=io, stderr=io))
        end
    end
end
println("Results: ", output)
