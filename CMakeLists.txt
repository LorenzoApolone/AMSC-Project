cmake_minimum_required(VERSION 3.14)
project(DPSO_Test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Find MPI ---
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

# --- DPSO library (the solver) ---
add_library(dpso STATIC
    src_dpso/dpso.cpp
    src_dpso/interfaces.cpp
    src_dpso/functions.cpp
)

target_include_directories(dpso PUBLIC
    ${CMAKE_SOURCE_DIR}/src_dpso
)

target_link_libraries(dpso PUBLIC
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
)

# --- Test executable ---
add_executable(test_dpso
    src_dpso/test_dpso.cpp
)

target_link_libraries(test_dpso PRIVATE
    dpso
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
)

# --- Optional: install target ---
install(TARGETS test_dpso DESTINATION bin)
