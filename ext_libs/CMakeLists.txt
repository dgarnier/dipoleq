include(FetchContent)

FetchContent_Declare(hdf4
                     GIT_REPOSITORY https://github.com/HDFGroup/hdf4.git
                     GIT_TAG        hdf-4_2_16
                    )

# set(HDF4_ENABLE_JPEG_LIB_SUPPORT OFF CACHE BOOL "" FORCE)
# set(BUILD_ZLIB_WITH_FETCHCONTENT ON CACHE BOOL "" FORCE)
set(HDF4_EXTERNALLY_CONFIGURED ON CACHE BOOL "" FORCE)
set(HDF4_EXPORTED_TARGETS "project-targets" CACHE STRING "" FORCE)
#set(BUILD_JPEG_WITH_FETCHCONTENT ON CACHE BOOL "" FORCE)
set(HDF4_BUILD_TOOLS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(hdf4)

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

#FetchContent_GetProperties(hdf4)


set(hdf4_INCLUDE_DIRS "${HDF4_MFHDF_LIBSRC_SOURCE_DIR}" "${HDF4_HDF_SRC_SOURCE_DIR}" "${HDF4_BINARY_DIR}" CACHE STRING "" FORCE)
set(hdf4_LIBRARIES "mfhdf-static" "hdf-static" CACHE STRING "" FORCE)

# hdf4_POPULATED, hdf4_SOURCE_DIR and hdf4_BINARY_DIR is defined by FetchContent
#find_package(hdf4)

#ExternalProject_Add(ClibPDF 
#                    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/clibpdf
#                    )

#ExternalProject_Get_Property(ClibPDF)

FetchContent_Declare(ClibPDF
                    OVERRIDE_FIND_PACKAGE
                    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/clibpdf
                    )

#ExternalProject_Get_Property(ClibPDF)
FetchContent_MakeAvailable(ClibPDF)
find_package(ClibPDF REQUIRED)

#dump_cmake_variables(ClibPDF)

#set(hdf4_INCLUDE_DIRS "${HDF4_MFHDF_LIBSRC_SOURCE_DIR}" "${HDF4_HDF_SRC_SOURCE_DIR}" "${HDF4_BINARY_DIR}" CACHE STRING "" FORCE)
#set(hdf4_LIBRARIES "mfhdf-static" "hdf-static" CACHE STRING "" FORCE)

# dump_cmake_variables()

