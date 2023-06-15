add_library(clibpdf 
    cpdfInit.c cpdfPreview.c cpdfAnnot.c cpdfUtil.c cpdfTime.c cpdfText.c
    cpdfRawPlot.c cpdfDomPlot.c cpdfDomain.c cpdfColor.c cpdfGstate.c
    cpdfMarker.c cpdfAxis.c cpdfFont.c cpdfTicInt.c cpdfMemBuf.c cpdfImage.c
    jpegsize.c cpdfArc.c cpdfImgIL.c cpdfOutline.c cpdfReadPFM.c cpdfReadPFB.c
    cpdfTextBox.c
)

add_library(clibpdf_zlib STATIC
    compress.c deflate.c adler32.c trees.c zutil.c
)

if (APPLE)
    target_compile_definitions(clibpdf PRIVATE MacOSX=1)
endif()

target_link_libraries(clibpdf PRIVATE clibpdf_zlib)
set_target_properties(clibpdf PROPERTIES PUBLIC_HEADER "cpdflib.h")
target_include_directories(clibpdf INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# dump_cmake_variables()