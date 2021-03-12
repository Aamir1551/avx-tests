cmake_minimum_required(VERSION 3.17)
project(simd_prac)

set(CMAKE_CXX_STANDARD 14)

add_executable(simd_prac main.cpp)
target_include_directories(simd_prac PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/glm")
target_compile_options(simd_prac PUBLIC  "-mavx" "-mavx2")
