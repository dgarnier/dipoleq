# dipoleq
cmake_minimum_required(VERSION 3.20)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(DipolEq VERSION 1.2 LANGUAGES C CXX)

include_directories("${CMAKE_SOURCE_DIR}/includes")
add_subdirectory(ext_libs EXCLUDE_FROM_ALL)
# for some reason, ext_libs can't export this info...
include(FindHDF4)

add_compile_definitions(DIPOLE=1 HDFOUTPUT=1)
add_subdirectory(fortran)
add_subdirectory(sources EXCLUDE_FROM_ALL)

add_executable(dipoleq)

target_sources(dipoleq
    PUBLIC
        SimDipEq.c)

target_compile_definitions(dipoleq PUBLIC PDFOUTPUT=1)

target_link_libraries(dipoleq PUBLIC DipolEqLib HDF4::HDF4 ${ClibPDF_LIBRARIES})
target_include_directories(dipoleq PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/includes")

install(TARGETS dipoleq DESTINATION bin)

