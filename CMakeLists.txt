# dipoleq
cmake_minimum_required(VERSION 3.24)

project(DipolEq VERSION 1.2 LANGUAGES C CXX)

include_directories("${CMAKE_SOURCE_DIR}/includes")
add_subdirectory(ext_libs)

add_compile_definitions(DIPOLE=1 HDFOUTPUT=1)
add_subdirectory(sources)

add_executable(dipoleq)

target_sources(dipoleq
    PUBLIC
        SimDipEq.c)

target_compile_definitions(dipoleq PUBLIC PDFOUTPUT=1)

target_link_libraries(dipoleq PUBLIC DipolEqLib ${hdf4_LIBRARIES})
target_include_directories(dipoleq PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/includes" "${hdf4_INCLUDE_DIRS}" "${HDF4_INSTALL_INCLUDE_DIR}")

