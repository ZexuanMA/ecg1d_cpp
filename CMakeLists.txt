cmake_minimum_required(VERSION 3.16)
project(ecg1d_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fetch Eigen
include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen)

# Main executable with all source files
add_executable(ecg1d
    main.cpp
    basis_params.cpp
    permutation.cpp
    pair_cache.cpp
    interaction_kernels.cpp
    hamiltonian.cpp
    derivatives.cpp
    observable_derivatives.cpp
    hamiltonian_gradient.cpp
    tdvp_solver.cpp
    svm.cpp
)

target_link_libraries(ecg1d PRIVATE Eigen3::Eigen)
target_compile_options(ecg1d PRIVATE -O2 -march=native)
