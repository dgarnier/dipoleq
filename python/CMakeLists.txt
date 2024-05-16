
Python_add_library(pydipoleq MODULE WITH_SOABI PyDipolEq.cpp)

target_compile_definitions(pydipoleq PRIVATE PYBIND11=1)
target_link_libraries(pydipoleq PRIVATE pydipoleqlib pybind11::module)
set_target_properties(pydipoleq PROPERTIES OUTPUT_NAME "_c")
