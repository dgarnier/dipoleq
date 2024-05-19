execute_process(COMMAND pybind11-config --cmakedir
    RESULT_VARIABLE PYBIND11_CONFIG_RESULT
    ERROR_QUIET
    OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT PYBIND11_CONFIG_RESULT EQUAL 0)
    message(STATUS "pybind11-config failed to execute")
else()
    message(STATUS "Found pybind11Config.cmake at ${PYBIND11_CMAKE_DIR}")
endif()
if(PYBIND11_CONFIG_RESULT EQUAL 0)
    set(pybind11_DIR ${PYBIND11_CMAKE_DIR} CACHE PATH "Path to pybind11Config.cmake")
endif()
