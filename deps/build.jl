# Pkg.build("DACE"): build the unmodified, pinned upstream C engine.
using Libdl
source = isempty(ARGS) ? joinpath(@__DIR__, "source-upstream") : abspath(only(ARGS))
revision = "368f6020342b0223ef1c05d50602bc54b7a8fe1f"
if !isdir(source)
    run(`git clone --no-checkout https://github.com/dacelib/dace $source`)
    run(`git -C $source checkout --detach $revision`)
end
readchomp(`git -C $source rev-parse HEAD`) == revision || error("Expected native revision $revision")
build = joinpath(@__DIR__, "build-upstream")
prefix = joinpath(@__DIR__, "usr")
if Sys.iswindows() && (Sys.which("gcc") === nothing || Sys.which("cmake") === nothing)
    include("build_windows.jl")
else
    generator = Sys.iswindows() ? ["-G", "MinGW Makefiles"] : String[]
    run(`cmake -S $(@__DIR__) -B $build $generator -DDACE_SOURCE=$source -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$prefix`)
    run(`cmake --build $build --target dacecore --parallel 4`)
    run(`cmake --install $build`)
end
println("Built ", joinpath(prefix, "lib", "libdacecore.$dlext"))
