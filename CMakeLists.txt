# CMake file for DipolEq - a dipole equilibrium solver
cmake_minimum_required(VERSION 3.26)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# dynamic versioning from git tags
include(get_cpm)
CPMAddPackage("gh:LecrisUT/CMakeExtraUtils@0.4.1")
list(APPEND CMAKE_MODULE_PATH "${CMakeExtraUtils_BINARY_DIR}")
include(DynamicVersion)
dynamic_version(PROJECT_PREFIX DipolEq_ 
                FALLBACK_VERSION 0.9.999
                ALLOW_FAILS
                )

project(DipolEq VERSION ${PROJECT_VERSION} LANGUAGES C CXX)
include(CTest)

option(BUILD_CLI "Build the CLI executable" ON)
option(BUILD_WITH_HDF5 "Build with HDF5" ON)
option(BUILD_WITH_HDF4 "Build with HDF4" OFF)
option(BUILD_HDF4_WITH_FETCHCONTENT "Build HDF4 with FetchContent" OFF)
option(BUILD_WITH_CLIBPDF "Build with C libPDF" ON)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if (NOT BUILD_CLI)
    set(BUILD_WITH_HDF5 OFF)
    set(BUILD_WITH_HDF4 OFF)
    set(BUILD_WITH_CLIBPDF OFF)
endif()

if (SKBUILD)
    set(BUILD_PYTHON_BINDINGS ON)
endif()

if (BUILD_PYTHON_BINDINGS)
    # maybe we should force make a virtualenv here
    # cmake 3.26 for changes to FindPython
    set(Python_FIND_VIRTUALENV "FIRST")
    set(Python_FIND_FRAMEWORK "LAST")
    find_package(Python GLOBAL REQUIRED 
        COMPONENTS 
            Interpreter Development.Module NumPy 
            ${SKBUILD_SABI_COMPONENT}
    )
endif()

add_subdirectory(ext_libs EXCLUDE_FROM_ALL)

include_directories("${CMAKE_SOURCE_DIR}/includes")

add_compile_definitions(DIPOLE=1)


add_subdirectory(sources EXCLUDE_FROM_ALL)

add_executable(dipoleq)
target_sources(dipoleq PUBLIC SimDipEq.c)
if (BUILD_WITH_CLIBPDF)
    target_compile_definitions(dipoleq PUBLIC PDFOUTPUT=1)
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
if (BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if (BUILD_CLI)
    add_subdirectory(Testing)
endif()

# for some reason.. TARGETS must install from top level
if (SKBUILD)
    message(STATUS "SKBUILD STATUS = ${SKBUILD_STATE}")
    install(TARGETS pydipoleq LIBRARY DESTINATION dipoleq COMPONENT python)
endif()

