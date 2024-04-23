include(FetchContent)
include(FindZLIB)
include(FindJPEG)

function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
        if (ARGV0)
            unset(MATCHED)
            string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
            if (NOT MATCHED)
                continue()
            endif()
        endif()
        message(STATUS "${_variableName}=${${_variableName}}")
    endforeach()
endfunction()

if(NOT JPEG_FOUND)
    include(ExternalProject)
    ExternalProject_Add(
        jpeg
        URL http://ijg.org/files/jpegsrc.v9e.tar.gz
        URL_HASH SHA256=4077d6a6a75aeb01884f708919d25934c93305e49f7e3f36db9129320e6f4f3d
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/jpeg-src"
        BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/jpeg-build"
        INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/jpeg-install"
        CONFIGURE_COMMAND "${CMAKE_CURRENT_BINARY_DIR}/jpeg-src/configure"
            "--prefix=${CMAKE_CURRENT_BINARY_DIR}/jpeg-install"
            "--enable-static"
            "--disable-shared"
            "--with-pic"
        BUILD_COMMAND "make"
        INSTALL_COMMAND "make" "install"
    )

    set(JPEG_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/jpeg-install/include" CACHE STRING " " FORCE)
    set(JPEG_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/jpeg-install/lib/libjpeg.a" CACHE STRING " " FORCE)

    add_library(LOCAL_JPEG STATIC IMPORTED GLOBAL)
    add_dependencies(LOCAL_JPEG jpeg)
    set_target_properties(LOCAL_JPEG PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${JPEG_INCLUDE_DIR}"
        IMPORTED_LOCATION "${JPEG_LIBRARY}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
    )
    add_library(JPEG::JPEG ALIAS LOCAL_JPEG)
    
endif()
# dump_cmake_variables(JPEG)

if (BUILD_WITH_HDF4)
    if(NOT BUILD_HDF4_WITH_FETCHCONTENT)
        find_package(HDF4)
        if(TARGET HDF4::HDF4)
            dump_cmake_variables(HDF4)
        else()
            set(BUILD_HDF4_WITH_FETCHCONTENT ON CACHE BOOL " " FORCE)
        endif()
    endif()

    if(BUILD_HDF4_WITH_FETCHCONTENT)
        # if you don't have it, build it
        FetchContent_Declare(HDF4
                            GIT_REPOSITORY https://github.com/HDFGroup/hdf4.git
                            GIT_TAG        hdf-4_2_16-2
                            )

        # set(HDF4_ENABLE_JPEG_LIB_SUPPORT OFF CACHE BOOL " " FORCE)
        # set(BUILD_ZLIB_WITH_FETCHCONTENT ON CACHE BOOL " " FORCE)
        set(HDF4_EXTERNALLY_CONFIGURED ON CACHE BOOL " " FORCE)
        set(HDF4_EXPORTED_TARGETS "project-targets" CACHE STRING " " FORCE)
        # set(BUILD_JPEG_WITH_FETCHCONTENT ON CACHE BOOL " " FORCE)
        set(HDF4_BUILD_TOOLS OFF CACHE BOOL " " FORCE)

        FetchContent_MakeAvailable(HDF4)

        set(HDF4_INCLUDE_DIRS "${HDF4_MFHDF_LIBSRC_SOURCE_DIR}" "${HDF4_HDF_SRC_SOURCE_DIR}" "${HDF4_BINARY_DIR}" CACHE STRING " " FORCE)
        set(HDF4_LIBRARIES "mfhdf-static" "hdf-static" CACHE STRING " " FORCE)

        add_library(hdf4 INTERFACE)
        set_target_properties(hdf4 PROPERTIES
                            INTERFACE_INCLUDE_DIRECTORIES "${HDF4_INCLUDE_DIRS}"
                            INTERFACE_LINK_LIBRARIES "${HDF4_LIBRARIES}"
                            )
        add_dependencies(hdf4 JPEG::JPEG)

        add_library(HDF4::HDF4 ALIAS hdf4)
        # dump_cmake_variables(HDF4)
    endif()
endif()

if (CMAKE_VERSION VERSION_LESS 3.24)
    FetchContent_Declare(ClibPDF
                        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/clibpdf
                        )
    FetchContent_MakeAvailable(ClibPDF)
else()
    FetchContent_Declare(ClibPDF
                        OVERRIDE_FIND_PACKAGE
                        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/clibpdf
                        )

    FetchContent_MakeAvailable(ClibPDF)
    find_package(ClibPDF REQUIRED)
endif()



