# Temporary source build until DACE_jll ships the evaluate binding.
# Run with the package project active; requires CMake, Ninja, and a C++ compiler.
using Pkg, TOML, Libdl
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
project = dirname(Base.active_project())
Pkg.instantiate(; allow_autoprecomp=false)
using libcxxwrap_julia_jll
cxxwrap = libcxxwrap_julia_jll.artifact_dir
Pkg.activate(; temp=true)
Pkg.add(PackageSpec(name="Eigen_jll", version="3.4.0"))
using Eigen_jll
eigen = Eigen_jll.artifact_dir
Pkg.activate(project)

root = joinpath(project, ".native")
source = joinpath(root, "source")
prefix = joinpath(root, "install")
revision = "bc26842ae5c6936f9e008a2dc0d9227a8d982acc"
if !isdir(source)
    mkpath(root)
    run(`git clone https://github.com/UoA-AstroGroup/dace.git $source`)
    run(`git -C $source checkout --detach $revision`)
end
@assert readchomp(`git -C $source rev-parse HEAD`) == revision
# Rename only the Julia bindings; the C++ evaluation API is unchanged.
wrapper = joinpath(source, "interfaces/julia/dace_julia.cxx")
write(wrapper, replace(read(wrapper, String), "mod.method(\"eval\"," => "mod.method(\"evaluate\","))
julia_prefix = dirname(Sys.BINDIR)
flags = Sys.iswindows() ? ["-DCMAKE_C_COMPILER=gcc", "-DCMAKE_CXX_COMPILER=g++"] : String[]
Sys.isapple() && push!(flags, "-DCMAKE_OSX_ARCHITECTURES=$(Sys.ARCH == :aarch64 ? "arm64" : "x86_64")")
run(`cmake -S $source -B $root/build -G Ninja -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=$prefix -DJulia_PREFIX=$julia_prefix
    -DJlCxx_DIR=$cxxwrap/lib/cmake/JlCxx -DEigen3_DIR=$eigen/share/eigen3/cmake
    -DWITH_JULIA=ON -DWITH_PTHREAD=ON -DWITH_ALGEBRAICMATRIX=ON
    -DWITH_EIGEN=ON -DCUSTOM_EXIT=ON $flags`)
run(`cmake --build $root/build --target install --parallel 2`)
library = joinpath(prefix, "lib", "libdace.$(Libdl.dlext)")
@assert isfile(library)
preferences_file = joinpath(project, "LocalPreferences.toml")
preferences = isfile(preferences_file) ? TOML.parsefile(preferences_file) : Dict{String,Any}()
get!(preferences, "DACE_jll", Dict{String,Any}())["libdace_path"] = library
open(io -> TOML.print(io, preferences), preferences_file, "w")
