cmake_minimum_required(VERSION 3.14)
project(GeneralizedKrawczykSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 自动下载Eigen库
include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW ON
    GIT_PROGRESS ON
)
FetchContent_MakeAvailable(eigen)

# 创建库
add_library(IntervalKrawczyk STATIC
    include/interval_krawczyk/KaucherInterval.h
    include/interval_krawczyk/IntervalVector.h
    include/interval_krawczyk/IntervalMatrix.h
    include/interval_krawczyk/GeneralizedKrawczykSolver.h
)

# 设置库的包含目录
target_include_directories(IntervalKrawczyk PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# 链接Eigen库
target_link_libraries(IntervalKrawczyk PUBLIC Eigen3::Eigen)

# 安装库
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

# 生成配置文件
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

# 添加测试
enable_testing()
add_subdirectory(tests)

# 添加演示程序
add_subdirectory(demos)
