add_library(DipolEqLib STATIC
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
    HDFOutput.c
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
    plot_contour.c PDFOutput.c
)

# TSPACK?
# target_compile_definitions(DipolEqLib PRIVATE _USE_TSPACK_=1)

target_link_libraries(DipolEqLib 
                        PRIVATE 
                            clibpdf
                            # tspack
                            HDF4::HDF4
                    )

target_include_directories(DipolEqLib 
                            PUBLIC 
                                ${CMAKE_SOURCE_DIR}/includes
                            )
