project(HFtest)

add_executable(${PROJECT_NAME} HFtest.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC
    gtest_main
    HFlib
)

add_test(
    NAME ${PROJECT_NAME}
    COMMAND ${PROJECT_NAME}
)