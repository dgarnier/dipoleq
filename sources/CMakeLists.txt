set(dipoleq_SOURCES)
list(APPEND dipoleq_SOURCES 
    AddCoilJ.c 
    AddShellJ.c 
    CDipoleIntStable.cc 
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
    GetFluxParameters.c
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
    list(APPEND dipoleq_LIBS hdf5::hdf5 hdf5::hdf5_hl)
    list(APPEND dipoleq_DEFS HDFOUTPUT=1)
endif()

if (BUILD_WITH_CLIBPDF)
    list(APPEND dipoleq_SOURCES plot_contour.c PDFOutput.c)
    list(APPEND dipoleq_LIBS clibpdf)
endif()

add_library(DipolEqLib STATIC ${dipoleq_SOURCES})

# TSPACK?
# target_compile_definitions(DipolEqLib PRIVATE _USE_TSPACK_=1)

target_link_libraries(DipolEqLib PRIVATE ${dipoleq_LIBS})

target_include_directories(DipolEqLib 
                            PUBLIC 
                                ${CMAKE_SOURCE_DIR}/includes
                            )
target_compile_definitions(DipolEqLib PRIVATE ${dipoleq_DEFS})

# Python bindings library for pydipoleq target
if (BUILD_PYTHON_BINDINGS)
    add_library(pydipoleqlib STATIC ${pydipoleqlib_SOURCES})
    target_include_directories(pydipoleqlib PUBLIC ${CMAKE_SOURCE_DIR}/includes)
    target_compile_definitions(pydipoleqlib PRIVATE PYBIND11=1)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_compile_options(pydipoleqlib PRIVATE -fPIC)
    endif()
endif()