using TOML, Printf
root = normpath(joinpath(@__DIR__, ".."))
folder = isempty(ARGS) ? joinpath(root, "results", "paired") : only(ARGS)
median(values) = (x = sort(values); (x[(length(x)+1)÷2] + x[(length(x)+2)÷2])/2)

for tag in ("113", "110")
    direct = [TOML.parsefile(joinpath(folder, "direct-$tag-$round.toml")) for round in 1:2]
    cxx = [TOML.parsefile(joinpath(folder, "cxxwrap-$tag-$round.toml")) for round in 1:2]
    println("Julia ", direct[1]["julia"])
    println("| Workload | CxxWrap (us) | Direct (us) | Direct/CxxWrap |")
    println("| --- | ---: | ---: | ---: |")
    for i in eachindex(direct[1]["rows"])
        drows = [run["rows"][i] for run in direct]
        crows = [run["rows"][i] for run in cxx]
        @assert all(row["name"] == drows[1]["name"] for row in vcat(drows, crows))
        d = median(vcat([r["seconds"] for r in drows]...))
        c = median(vcat([r["seconds"] for r in crows]...))
        @printf("| %s | %.3f | %.3f | %.3f |\n", drows[1]["name"], 1e6*c, 1e6*d, d/c)
    end
    println()
end
