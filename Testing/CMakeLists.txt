include(CTest)

add_test(NAME run_dipoleq COMMAND dipoleq -f ${CMAKE_CURRENT_SOURCE_DIR}/beta1m.in)
add_test(NAME run_dipoleq_7 
         COMMAND dipoleq -f ${CMAKE_CURRENT_SOURCE_DIR}/../python/tests/beta1_psinpeak.in
        )
add_test(NAME create_gfile COMMAND python ${dipoleq_geqdsk_SCRIPT} beta1.h5)
set_tests_properties(create_gfile PROPERTIES DEPENDS run_dipoleq)