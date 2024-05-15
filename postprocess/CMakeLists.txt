# this works, even though its not at the top level...
# I really don't understand cmake...

install(PROGRAMS h5togeqdsk.py 
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME "dipoleq_geqdsk"
    COMPONENT cli
)

set(dipoleq_geqdsk_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/h5togeqdsk.py PARENT_SCOPE)