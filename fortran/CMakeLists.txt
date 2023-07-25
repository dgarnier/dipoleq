set(CMAKE_C_STANDARD 99)

add_library(tspack STATIC
    tspack.c
    stcom_com.c
    tspack.h
)

set_target_properties(tspack PROPERTIES PUBLIC_HEADER "tspack.h")
target_include_directories(tspack INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})