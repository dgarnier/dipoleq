# dyanmic versioning, in c/c++ files
# see https://github.com/LecrisUT/CMakeExtraUtils/blob/main/cmake/DynamicVersion.md
configure_file(../includes/version.hpp.in version.hpp)
configure_file(../includes/version.h.in version.h)
add_library(version_lib INTERFACE version.hpp version.h)
target_include_directories(version_lib INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
add_dependencies(version_lib DipolEq_Version)
set_property(SOURCE version.hpp.in APPEND PROPERTY
		OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/.version)
set_property(SOURCE version.hpp.in APPEND PROPERTY
		OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/.version)


set(dipoleq_SOURCES)
list(APPEND dipoleq_SOURCES
    AddCoilJ.c
    AddShellJ.c
    CDipoleIntStable.cc
    CDipoleStablePsiN.cc
    CDipoleStd.cc
    CPlasmaModel.cc
    DelChiSqr.c
    FileInput.cc
    FileOutput.cc
    FindJ.cc
    FindMeasFit.c
    Find_ShellCurrent.c
    Find_dJdy.c
    GetFluxMoments.c
    GetFluxParameters.cc
    GetPlasmaParameters.cc
    InitJ.c
    J_DipoleStd.c
    J_IsoNoFlow.c
    J_Std.c
    LeastSquares.c
    LoadBndryGreens.c
    LoadMeasGreens.c
    LoadShellGreens.c
    PlasmaBoundary.c
    PsiBoundary.c
    Restart.c
    SVDFit.c
    bcucof.c
    bcuint.c
    coil.c
    contour.c
    dUnkn.c
    fpoly.c
    green.c
    interpolate.c
    limiter.c
    ludcmp.c
    meas_J0.c
    meas_bp.c
    meas_bpangle.c
    meas_circle.c
    meas_coilcur.c
    meas_flux.c
    meas_mag_Xxx.c
    meas_plasmacur.c
    meas_pnorm.c
    meas_ppsix.c
    meas_press.c
    meas_saddle.c
    measurement.c
    multitask.c
    nrRomberg.c
    nrSpline.c
    nrutil.c
    plasma.cc
    polygon.c
    psigrid.c
    rolldown.c
    separatrix.c
    shell.c
    stdio_dmatrix.c
    tokamak.c
    tokgreen.c
)

set(pydipoleqlib_SOURCES ${dipoleq_SOURCES})

set(dipoleq_LIBS)
set(dipoleq_DEFS)

if (BUILD_WITH_HDF4)
    list(APPEND dipoleq_SOURCES HDFOutput.c)
    list(APPEND dipoleq_LIBS HDF4::HDF4)
    list(APPEND dipoleq_DEFS HDFOUTPUT=1)
endif()

if (BUILD_WITH_HDF5)
    list(APPEND dipoleq_SOURCES HDF5Output.c)

    if (TARGET HDF5::HDF5)
        message(STATUS "HDF5::HDF5 found")
        list(APPEND dipoleq_LIBS HDF5::HDF5)
    else()
        if (TARGET hdf5::hdf5)
            message(STATUS "hdf5::hdf5 found")
            list(APPEND dipoleq_LIBS hdf5::hdf5)
        else()
            message(STATUS "HDF5::HDF5 not found")
            message(STATUS "hdf5::hdf5 not found")
        endif()
    endif()

    if (TARGET hdf5::hdf5_hl)
        message(STATUS "hdf5::hdf5_hl found")
        list(APPEND dipoleq_LIBS hdf5::hdf5_hl)
    else()
        message(STATUS "HDF5::HL not found")
    endif()
    list(APPEND dipoleq_DEFS HDFOUTPUT=1)
endif()

if (BUILD_WITH_CLIBPDF)
    list(APPEND dipoleq_SOURCES plot_contour.c PDFOutput.c)
    list(APPEND dipoleq_LIBS clibpdf)
endif()

add_library(DipolEqLib STATIC EXCLUDE_FROM_ALL ${dipoleq_SOURCES})

if (LINUX AND (DEBUG_MODE OR SANITIZER_ENABLED))
    target_link_libraries(DipolEqLib PUBLIC --fsanitize=address)
    target_compile_definitions(DipolEqLib PUBLIC --fsanaitize=address)
endif()

# TSPACK?
# target_compile_definitions(DipolEqLib PRIVATE _USE_TSPACK_=1)

target_link_libraries(DipolEqLib PRIVATE ${dipoleq_LIBS}
                                 PUBLIC version_lib)

target_include_directories(DipolEqLib
                            PUBLIC
                                ${CMAKE_SOURCE_DIR}/includes
                            )
target_compile_definitions(DipolEqLib PRIVATE ${dipoleq_DEFS})

# Python bindings library for pydipoleq target
if (BUILD_PYTHON_BINDINGS)
    add_library(pydipoleqlib EXCLUDE_FROM_ALL STATIC ${pydipoleqlib_SOURCES})
    target_include_directories(pydipoleqlib PUBLIC ${CMAKE_SOURCE_DIR}/includes)
    target_compile_definitions(pydipoleqlib PRIVATE PYBIND11=1)
    target_link_libraries(pydipoleqlib PUBLIC version_lib)
    set_target_properties(pydipoleqlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
