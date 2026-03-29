
# python is already installed at this point
# and we have pybind11 and/or nanobind installed

if(BUILD_USE_NANOBIND)

    find_package(nanobind CONFIG REQUIRED)
    nanobind_add_module(pydipoleq NB_STATIC STABLE_ABI
        NbDipolEq.cpp
    )
    target_compile_definitions(pydipoleq PRIVATE NANOBIND=1)
    target_link_libraries(pydipoleq PRIVATE pydipoleqlib)
    set_target_properties(pydipoleq PROPERTIES
        OUTPUT_NAME "core"
        POSITION_INDEPENDENT_CODE ON
    )
    # nanobind has cmake stub generation command
    nanobind_add_stub(
        pydipoleq_stubs
        MODULE core
        OUTPUT core.pyi
        PYTHON_PATH $<TARGET_FILE_DIR:pydipoleq>
        DEPENDS pydipoleq
    )

    if(NOT WIN32)
        add_custom_command(TARGET pydipoleq_stubs POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
                $<TARGET_FILE:pydipoleq> $<TARGET_FILE_DIR:pydipoleq>/core.pyi
                ${CMAKE_CURRENT_SOURCE_DIR}/dipoleq
            COMMENT "Copying nanobind module to dipoleq\n"
        )
    endif()

    install(TARGETS pydipoleq
        LIBRARY DESTINATION DipolEq)

else()
    # ── pybind11 binding ─────────────────────────────────────────────────────
    include(FindPyBind11)
    find_package(pybind11 REQUIRED)
    if(NOT "${SKBUILD_SABI_COMPONENT}" STREQUAL "")
        python_add_library(pydipoleq MODULE WITH_SOABI USE_SABI 3.10 PyDipolEq.cpp)
    else()
        python_add_library(pydipoleq MODULE WITH_SOABI PyDipolEq.cpp)
    endif()
    target_compile_definitions(pydipoleq PRIVATE PYBIND11=1)
    target_link_libraries(pydipoleq PRIVATE pydipoleqlib pybind11::module)
    set_target_properties(pydipoleq PROPERTIES
        OUTPUT_NAME "core"
        POSITION_INDEPENDENT_CODE ON
    )

    if(NOT WIN32)
        add_custom_command(TARGET pydipoleq POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
                $<TARGET_FILE:pydipoleq>
                ${CMAKE_CURRENT_SOURCE_DIR}/dipoleq
            COMMENT "Copying pybind11 module to dipoleq\n"
        )
        add_custom_command(TARGET pydipoleq POST_BUILD
            COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/dev/stubgen.py
                    $<TARGET_FILE:pydipoleq>
                    ${CMAKE_CURRENT_SOURCE_DIR}/dipoleq
            COMMENT "Generating Python stubs\n"
        )
    endif()
endif()


# this works, even though its not at the top level...
# I really don't understand cmake...

install(PROGRAMS dipoleq/h5togeqdsk.py
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME "dipoleq_geqdsk"
    COMPONENT cli
)

set(dipoleq_geqdsk_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/dipoleq/h5togeqdsk.py PARENT_SCOPE)
