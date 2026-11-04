cmake_minimum_required(VERSION 3.20)
project(kpm_kubo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

enable_testing()

# Eigen backs the ED Lehmann oracle only; the KPM path (R2+) stays
# dependency-free, same idiom as the sibling repos.
find_package(Eigen3 REQUIRED NO_MODULE)

# DESIGN.md section 6 decision, taken at R1: kpm-tightbinding is a build
# dependency at a pinned commit. Its own CMake project requires CUDA
# (project(... CXX CUDA) and find_package(CUDAToolkit REQUIRED)), which a
# host-only module cannot take as a configure requirement, so the pinned
# sources are populated without add_subdirectory and the host core is
# compiled directly. Nothing is copied into this tree; the clone lives in
# the build directory; that repo is never modified.
# Default source is the public GitHub repo; for local development point
# KPM_TIGHTBINDING_REPO at a local checkout (any git source FetchContent
# accepts), e.g. -DKPM_TIGHTBINDING_REPO=../kpm-tightbinding.
set(KPM_TIGHTBINDING_REPO "https://github.com/FireHorseAutumn/kpm-tightbinding.git"
    CACHE STRING "kpm-tightbinding git URL or local checkout to build against")
set(KPM_TIGHTBINDING_PIN "c596c86a0e343be25e735a8be2433ee163ac0809")
include(FetchContent)
FetchContent_Declare(kpm_tightbinding
  GIT_REPOSITORY "${KPM_TIGHTBINDING_REPO}"
  GIT_TAG "${KPM_TIGHTBINDING_PIN}"
  SOURCE_SUBDIR __populate_only__)
FetchContent_MakeAvailable(kpm_tightbinding)

add_library(kpm_tightbinding_core STATIC
  ${kpm_tightbinding_SOURCE_DIR}/src/lattice.cpp)
target_include_directories(kpm_tightbinding_core PUBLIC
  ${kpm_tightbinding_SOURCE_DIR}/include)
target_compile_options(kpm_tightbinding_core PRIVATE -Wall -Wextra)

add_library(kubo STATIC src/kubo.cpp)
target_include_directories(kubo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(kubo PUBLIC kpm_tightbinding_core PRIVATE Eigen3::Eigen)
target_compile_options(kubo PRIVATE -Wall -Wextra)

# A CUDA target exists only if the R6 GPU decision is made.
add_executable(test_stub tests/test_stub.cpp)
target_link_libraries(test_stub PRIVATE kubo)
target_compile_options(test_stub PRIVATE -Wall -Wextra)
add_test(NAME stub COMMAND test_stub)

add_executable(test_r1 tests/test_r1.cpp)
target_link_libraries(test_r1 PRIVATE kubo Eigen3::Eigen)
target_compile_options(test_r1 PRIVATE -Wall -Wextra)
target_compile_definitions(test_r1 PRIVATE
  KPM_TB_DATA_DIR="${kpm_tightbinding_SOURCE_DIR}/data")
add_test(NAME r1 COMMAND test_r1)
