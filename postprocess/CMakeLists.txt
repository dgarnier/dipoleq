include(GNUInstallDirs)

install(PROGRAMS dipoleq_h5togeqdsk.py 
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME "dipoleq_geqdsk")

set(dipoleq_geqdsk_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/dipoleq_h5togeqdsk.py PARENT_SCOPE)