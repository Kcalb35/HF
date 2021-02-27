cmake_minimum_required(VERSION 3.10)

project(HF)

add_executable(${PROJECT_NAME} main.cpp)


enable_testing()

add_subdirectory(test)
add_subdirectory(lib/googletest)
add_subdirectory(src)

target_link_libraries(${PROJECT_NAME} HFlib)