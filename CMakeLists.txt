# CMake file for DipolEq - a dipole equilibrium solver
cmake_minimum_required(VERSION 3.20)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(DipolEq VERSION 1.2 LANGUAGES C CXX)
include(CTest)

option(BUILD_WITH_HDF5 "Build with HDF5" ON)
option(BUILD_WITH_HDF4 "Build with HDF4" OFF)
option(BUILD_HDF4_WITH_FETCHCONTENT "Build HDF4 with FetchContent" OFF)
option(BUILD_WITH_CLIBPDF "Build with C libPDF" ON)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

add_subdirectory(ext_libs EXCLUDE_FROM_ALL)

include_directories("${CMAKE_SOURCE_DIR}/includes")

add_compile_definitions(DIPOLE=1)
if (BUILD_WITH_HDF4)
    include(FindHDF4)
endif()
if (BUILD_WITH_HDF5)
    find_package(HDF5 COMPONENTS C HL)
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options(-fsanitize=address)
endif()

# add_subdirectory(fortran)
add_subdirectory(sources EXCLUDE_FROM_ALL)

add_executable(dipoleq)
target_sources(dipoleq PUBLIC SimDipEq.c)
if (BUILD_WITH_CLIBPDF)
    target_compile_definitions(dipoleq PUBLIC PDFOUTPUT=1)
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_options(dipoleq PRIVATE -fsanitize=address)
    target_link_options(dipoleq PRIVATE -fsanitize=address)
endif()

target_link_libraries(dipoleq PUBLIC DipolEqLib)
target_include_directories(dipoleq PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/includes")

if (SKBUILD)
    set(CMAKE_INSTALL_BINDIR ${SKBUILD_SCRIPTS_DIR})
else()
    set(CMAKE_INSTALL_BINDIR bin)
endif()

install(TARGETS dipoleq DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT cli)

add_subdirectory(postprocess)
add_subdirectory(python)
add_subdirectory(Testing)

# for some reason.. TARGETS must install from top level
if (SKBUILD)
    message(STATUS "SKBUILD STATUS = ${SKBUILD_STATE}")
    install(TARGETS pydipoleq LIBRARY DESTINATION dipoleq COMPONENT python)
endif()



