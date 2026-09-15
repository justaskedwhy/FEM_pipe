cmake_minimum_required(VERSION 3.16)
project(fem_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(pybind11 2.11 REQUIRED)

pybind11_add_module(femcore
    cpp/bindings/bindings.cpp
    cpp/src/solver.cpp
    cpp/src/elements_t3.cpp
    cpp/src/elements_q4.cpp
    cpp/src/elements_t4.cpp
    cpp/src/elements_h8.cpp
)

target_include_directories(femcore PRIVATE cpp/include)
target_link_libraries(femcore PRIVATE Eigen3::Eigen)
set_target_properties(femcore PROPERTIES INTERPROCEDURAL_OPTIMIZATION OFF)

install(TARGETS femcore DESTINATION .)