cmake_minimum_required(VERSION 3.10)

project(blonde_mapper VERSION 1.0.0)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/version.h
)

add_subdirectory(external/bioparser)

file(GLOB SRC_FILES src/*.cpp)
add_executable(blonde_mapper ${SRC_FILES})

target_include_directories(blonde_mapper PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    include
    external/bioparser/include
)

target_link_libraries(blonde_mapper PRIVATE bioparser)
