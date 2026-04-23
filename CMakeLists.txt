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
    src/basis_params.cpp
    src/permutation.cpp
    src/pair_cache.cpp
    src/interaction_kernels.cpp
    src/hamiltonian.cpp
    src/derivatives.cpp
    src/observable_derivatives.cpp
    src/hamiltonian_gradient.cpp
    src/tdvp_solver.cpp
    src/svm.cpp
    src/kicked_evolution.cpp
    src/kicked_exact.cpp
    src/kick_operator.cpp
    src/observables.cpp
    src/target_fitting.cpp
    src/realtime_tdvp.cpp
    src/h_squared.cpp
    src/rothe_operator.cpp
    src/rothe_step.cpp
)

target_include_directories(ecg1d PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(ecg1d PRIVATE Eigen3::Eigen)
target_compile_options(ecg1d PRIVATE -O2 -march=native)

# Step B test: orthogonalized subspace prototype
add_executable(stepb_test
    stepb_test.cpp
    src/basis_params.cpp
    src/permutation.cpp
    src/pair_cache.cpp
    src/interaction_kernels.cpp
    src/hamiltonian.cpp
    src/derivatives.cpp
    src/observable_derivatives.cpp
    src/hamiltonian_gradient.cpp
    src/tdvp_solver.cpp
    src/svm.cpp
    src/kicked_evolution.cpp
    src/kicked_exact.cpp
    src/kick_operator.cpp
    src/observables.cpp
    src/target_fitting.cpp
    src/realtime_tdvp.cpp
    src/h_squared.cpp
    src/rothe_operator.cpp
    src/rothe_step.cpp
)
target_include_directories(stepb_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(stepb_test PRIVATE Eigen3::Eigen)
target_compile_options(stepb_test PRIVATE -O2 -march=native)

# Simplified learning executable: only imaginary-time harmonic + grid kicked rotor
add_executable(ecg1d_simple
    main_simple.cpp
    src/basis_params.cpp
    src/permutation.cpp
    src/pair_cache.cpp
    src/interaction_kernels.cpp
    src/hamiltonian.cpp
    src/derivatives.cpp
    src/observable_derivatives.cpp
    src/hamiltonian_gradient.cpp
    src/tdvp_solver.cpp
    src/svm.cpp
    src/kick_operator.cpp
    src/target_fitting.cpp
    src/realtime_tdvp.cpp
    src/h_squared.cpp
    src/rothe_operator.cpp
    src/rothe_step.cpp
)
target_include_directories(ecg1d_simple PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(ecg1d_simple PRIVATE Eigen3::Eigen)
target_compile_options(ecg1d_simple PRIVATE -O2 -march=native)
