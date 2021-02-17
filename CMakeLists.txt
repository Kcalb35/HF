cmake_minimum_required(VERSION 3.10)

project(HF)

add_executable(${PROJECT_NAME} main.cpp)

set(SOURCES
    src/basis.cpp
    src/input.cpp
    src/Intergal.cpp
)

add_library(HFlib STATIC ${SOURCES})

target_link_libraries(${PROJECT_NAME} HFlib)

enable_testing()

add_subdirectory(test)
add_subdirectory(lib/gtest)