cmake_minimum_required(VERSION 3.14)
project(GeneralizedKrawczykSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Auto-download dependencies
include(FetchContent)

# Eigen library
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW ON
    GIT_PROGRESS ON
)

# Google Test library
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
    GIT_SHALLOW ON
    GIT_PROGRESS ON
)

# Make dependencies available
FetchContent_MakeAvailable(eigen googletest)

# Create library
add_library(IntervalKrawczyk STATIC
    include/interval_krawczyk/KaucherInterval.h
    include/interval_krawczyk/IntervalVector.h
    include/interval_krawczyk/IntervalMatrix.h
    include/interval_krawczyk/GeneralizedKrawczykSolver.h
)

# Set library include directories
target_include_directories(IntervalKrawczyk PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Link Eigen library
target_link_libraries(IntervalKrawczyk PUBLIC Eigen3::Eigen)

# Install library
install(TARGETS IntervalKrawczyk
    EXPORT IntervalKrawczykTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(EXPORT IntervalKrawczykTargets
    FILE IntervalKrawczykTargets.cmake
    NAMESPACE IntervalKrawczyk::
    DESTINATION lib/cmake/IntervalKrawczyk
)

install(DIRECTORY include/interval_krawczyk
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
)

# Generate configuration files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/IntervalKrawczykConfigVersion.cmake"
    VERSION 1.0.0
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/IntervalKrawczykConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/IntervalKrawczykConfig.cmake"
    INSTALL_DESTINATION lib/cmake/IntervalKrawczyk
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/IntervalKrawczykConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/IntervalKrawczykConfigVersion.cmake"
    DESTINATION lib/cmake/IntervalKrawczyk
)

# Add tests
enable_testing()
add_subdirectory(tests)

# Add demos
add_subdirectory(demos)
