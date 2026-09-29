# Windows fallback for deps/build.jl: WSL Ubuntu with MinGW-w64 POSIX compilers.
wslpath(path) = readchomp(`wsl -d Ubuntu -- wslpath -a -u $(replace(abspath(path), '\\' => '/'))`)
cmake = wslpath(@__DIR__)
build_wsl = wslpath(joinpath(@__DIR__, "build-upstream-windows"))
prefix_wsl = wslpath(prefix)
native = wslpath(source)
toolchain = wslpath(joinpath(@__DIR__, "mingw-toolchain.cmake"))
run(`wsl -d Ubuntu -- cmake -S $cmake -B $build_wsl -DDACE_SOURCE=$native -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=$prefix_wsl -DCMAKE_TOOLCHAIN_FILE=$toolchain`)
run(`wsl -d Ubuntu -- cmake --build $build_wsl --target dacecore --parallel 4`)
run(`wsl -d Ubuntu -- cmake --install $build_wsl`)
