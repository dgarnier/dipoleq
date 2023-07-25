include(FetchContent)
include(FindZLIB)
include(FindJPEG)
include(FindHDF4)

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

#find_package(HDF4)
if(TARGET HDF4::HDF4)
    dump_cmake_variables(HDF4)
#    add_library(HDF4::HDF4 ALIAS HDF4_TARGETS)
else()
#if(NOT HDF4_FOUND)
    # if you don't have it, build it
    FetchContent_Declare(HDF4
                        GIT_REPOSITORY https://github.com/HDFGroup/hdf4.git
                        GIT_TAG        hdf-4_2_16-2
                        )

    # set(HDF4_ENABLE_JPEG_LIB_SUPPORT OFF CACHE BOOL "" FORCE)
    # set(BUILD_ZLIB_WITH_FETCHCONTENT ON CACHE BOOL "" FORCE)
    set(HDF4_EXTERNALLY_CONFIGURED ON CACHE BOOL "" FORCE)
    set(HDF4_EXPORTED_TARGETS "project-targets" CACHE STRING "" FORCE)
    # set(BUILD_JPEG_WITH_FETCHCONTENT ON CACHE BOOL "" FORCE)
    set(HDF4_BUILD_TOOLS OFF CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(HDF4)

    set(HDF4_INCLUDE_DIRS "${HDF4_MFHDF_LIBSRC_SOURCE_DIR}" "${HDF4_HDF_SRC_SOURCE_DIR}" "${HDF4_BINARY_DIR}" CACHE STRING "" FORCE)
    set(HDF4_LIBRARIES "mfhdf-static" "hdf-static" CACHE STRING "" FORCE)

    add_library(hdf4 INTERFACE)
    set_target_properties(hdf4 PROPERTIES
                        INTERFACE_INCLUDE_DIRECTORIES "${HDF4_INCLUDE_DIRS}"
                        INTERFACE_LINK_LIBRARIES "${HDF4_LIBRARIES}"
                        )

    add_library(HDF4::HDF4 ALIAS hdf4)
    # dump_cmake_variables(HDF4)
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

# dump_cmake_variables()

