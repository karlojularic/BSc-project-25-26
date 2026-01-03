cmake_minimum_required(VERSION 3.10)

project(blonde_mapper VERSION 1.0.0)

configure_file(
    version.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/version.h
)

add_subdirectory(external/bioparser)

set(SRC_FILES
    src/main.cpp
    src/team_alignment.cpp
    src/team_common.cpp
    src/team_minimizers.cpp # DODANO: izvorna datoteka za minimizere
)

add_executable(blonde_mapper ${SRC_FILES})

target_include_directories(blonde_mapper PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    include
)

target_link_libraries(blonde_mapper PRIVATE bioparser)

enable_testing()

find_package(GTest REQUIRED)

# DODANO: test/minimizer_tests.cpp u popis izvornih datoteka za testove
add_executable(blonde_mapper_tests 
    test/alignment_tests.cpp 
    test/minimizer_tests.cpp
)

target_include_directories(blonde_mapper_tests PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    include
    src
)

target_sources(blonde_mapper_tests PRIVATE
    src/team_alignment.cpp
    src/team_common.cpp
    src/team_minimizers.cpp # DODANO: izvorna datoteka za minimizere u testovima
)

target_link_libraries(blonde_mapper_tests PRIVATE
    GTest::gtest
    GTest::gtest_main
    bioparser
)

add_test(NAME blonde_mapper_tests COMMAND blonde_mapper_tests)