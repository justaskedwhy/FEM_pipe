cmake_minimum_required(VERSION 3.16)
project(fem_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Eigen3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    include(FetchContent)
    set(EIGEN_VERSION_TAG 3.4.0)
    message(STATUS "Eigen3 not found on system; downloading ${EIGEN_VERSION_TAG} via FetchContent")
    FetchContent_Declare(eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION_TAG}/eigen-${EIGEN_VERSION_TAG}.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen)
endif()
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