cmake_minimum_required(VERSION 4.0)

project(blonde_mapper VERSION 1.0.0)

add_executable(blonde_mapper blonde_mapper.cpp)

target_include_directories(blonde_mapper PRIVATE ${CMAKE_CURRENT_BINARY_DIR})


configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/version.h
)