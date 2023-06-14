# dipoleq
cmake_minimum_required(VERSION 3.24)

project(DipolEq VERSION 1.2 LANGUAGES C CXX)

include_directories("${CMAKE_SOURCE_DIR}/includes")
add_compile_definitions(DIPOLE=1)
add_subdirectory(sources)

add_executable(dipoleq)

target_sources(dipoleq
    PUBLIC
        SimDipEq.c)
target_link_libraries(dipoleq PUBLIC DipolEqLib)
target_include_directories(dipoleq PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/includes")