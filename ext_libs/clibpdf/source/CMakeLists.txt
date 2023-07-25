add_library(clibpdf STATIC
    cpdfInit.c cpdfPreview.c cpdfAnnot.c cpdfUtil.c cpdfTime.c cpdfText.c
    cpdfRawPlot.c cpdfDomPlot.c cpdfDomain.c cpdfColor.c cpdfGstate.c
    cpdfMarker.c cpdfAxis.c cpdfFont.c cpdfTicInt.c cpdfMemBuf.c cpdfImage.c
    jpegsize.c cpdfArc.c cpdfImgIL.c cpdfOutline.c cpdfReadPFM.c cpdfReadPFB.c
    cpdfTextBox.c
)

# don't use old broken ZLIB
#add_library(clibpdf_zlib STATIC
#    compress.c deflate.c adler32.c trees.c zutil.c
#)
#target_link_libraries(clibpdf PRIVATE clibpdf_zlib)

include(FindZLIB)

if (APPLE)
    target_compile_definitions(clibpdf PRIVATE MacOSX=1)
elseif(CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_compile_definitions(clibpdf PRIVATE Linux=1)
endif()

target_link_libraries(clibpdf PRIVATE ZLIB::ZLIB)
set_target_properties(clibpdf PROPERTIES PUBLIC_HEADER "cpdflib.h")
target_include_directories(clibpdf INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# dump_cmake_variables()