project(HFlib)

find_package(GSL REQUIRED)

set(SOURCES
    basis.cpp
    input.cpp
    Integral.cpp
)

add_library(${PROJECT_NAME} STATIC ${SOURCES})

target_link_libraries(${PROJECT_NAME} GSL::gsl)