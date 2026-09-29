cmake_minimum_required(VERSION 3.13)
project(DACEDirectCore LANGUAGES C CXX)
set(WITH_PTHREAD OFF CACHE BOOL "" FORCE)
set(DACE_MEMORY_MODEL DYNAMIC CACHE STRING "" FORCE)
# MinGW declares optional Annex K functions that its CRT does not provide.
# Keep this platform configuration in our build recipe, not an upstream fork.
if(MINGW)
  set(HAVE_SAFE_STRINGS FALSE CACHE INTERNAL "" FORCE)
endif()
add_subdirectory("${DACE_SOURCE}" dace EXCLUDE_FROM_ALL)

# Build only DACE's C engine; none of the C++ interface targets are linked.
add_library(dacecore SHARED abi.c $<TARGET_OBJECTS:c_core> $<TARGET_OBJECTS:c_core_contrib>)
target_include_directories(dacecore PRIVATE "${DACE_SOURCE}/core/include" "${CMAKE_CURRENT_BINARY_DIR}/dace/core/include")
if(WIN32)
  target_compile_definitions(dacecore PRIVATE "DACE_API=__declspec(dllexport)")
else()
  target_link_libraries(dacecore PRIVATE m)
endif()
install(TARGETS dacecore LIBRARY DESTINATION lib RUNTIME DESTINATION lib)
install(FILES "${DACE_SOURCE}/LICENSE" "${DACE_SOURCE}/NOTICE" DESTINATION share/licenses/DACE)
