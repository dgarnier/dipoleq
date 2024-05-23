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

include(FetchContent)
if(BUILD_WITH_HDF4 OR BUILD_WITH_HDF5)
    include(FindZLIB)

    if(NOT ZLIB_FOUND)
        include(get_cpm)
        CPMAddPackage(
            NAME zlib-cmake
            GITHUB_REPOSITORY jimmy-park/zlib-cmake
            GIT_TAG 1.3.1
            OPTIONS "ZLIB_USE_STATIC_LIBS"
        )
    endif()
endif()

if(BUILD_WITH_CLIBPDF OR BUILD_WITH_HDF5)
    include(FindJPEG)

    if (NOT JPEG_FOUND)
        include(get_cpm)
        CPMAddPackage(
            NAME jpeg-turbo
            GITHUB_REPOSITORY libjpeg-turbo/libjpeg-turbo
            GIT_TAG 3.0.3
            OPTIONS "ENABLE_SHARED=OFF"
        )
        add_library(JPEG::JPEG ALIAS jpeg-static)
    endif()

    if (NOT TARGET JPEG::JPEG)
        # not likely to have a CPMAddPackage for jpeg
        # as it doesn't use cmake
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
endif()

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

if (BUILD_WITH_HDF5)
    if (NOT LINUX)
        # static is what we want, but its broken on linux
        set(HDF5_USE_STATIC_LIBRARIES ON CACHE BOOL " " FORCE)
    endif()
    find_package(HDF5 GLOBAL COMPONENTS C HL)
    if(NOT HDF5_FOUND)
        include(get_cpm)
        CPMAddPackage(
            NAME hdf5
            GITHUB_REPOSITORY HDFGroup/hdf5
            GIT_TAG hdf5_1.14.4
            OPTIONS
                "BUILD_SHARED_LIBS OFF"
                "HDF5_BUILD_TOOLS OFF"
                "HDF5_BUILD_EXAMPLES OFF"
                "HDF5_BUILD_HL_LIB ON"
                "HDF5_BUILD_CPP_LIB OFF"
                "HDF5_BUILD_FORTRAN OFF"
                "HDF5_BUILD_JAVA OFF"
                "HDF5_BUILD_TESTING OFF"
                "HDF5_ENABLE_PARALLEL OFF"
                "HDF5_ENABLE_Z_LIB_SUPPORT ON"
                "HDF5_ENABLE_SZIP_SUPPORT OFF"
                "HDF5_ENABLE_FILTERS OFF"
                "HDF5_ENABLE_DEPRECATED_SYMBOLS OFF"
                "HDF5_ENABLE_CLEAR_FILE_BUFFER ON"
                "HDF5_ENABLE_LOGGING ON"
                "HDF5_ENABLE_USING_MEMCHECKER ON"
                "HDF5_ENABLE_USING_CODE_COVERAGE ON"
                "HDF5_ENABLE_USING_FI ON"
                "HDF5_EXTERNALLY_CONFIGURED ON"
        )

        # yes, this looks weird, but it matches
        # what needs happening from FindHDF5.cmake
        add_library(HDF5::HDF5 ALIAS hdf5-static)
        add_library(hdf5::hdf5_hl ALIAS hdf5_hl-static)

    endif()

endif()

if (BUILD_PYTHON_BINDINGS)
    # Find the module development requirements (requires FindPython from 3.17 or
    # scikit-build-core's built-in backport)

    # just look to install certain components.. python itself is 
    # previously found.

    include(FindPyBind11)
    find_package(pybind11 GLOBAL CONFIG QUIET)
    if(NOT pybind11_FOUND)
        include(get_cpm)
        CPMAddPackage(
            NAME pybind11
            VERSION 2.12.0
            GITHUB_REPOSITORY pybind/pybind11
            OPTIONS
                "PYBIND11_INSTALL ON"
                "PYBIND11_TEST OFF"
                "PYBIND11_FINDPYTHON OFF"
        )
        find_package(pybind11 GLOBAL REQUIRED)
    endif()
endif()

if(BUILD_WITH_CLIBPDF)
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
endif()