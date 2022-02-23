cmake_minimum_required(VERSION 3.12)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_C_COMPILER clang++-12)
set(CMAKE_CXX_COMPILER clang++-12)

project(spiral VERSION 0.1 LANGUAGES CXX)

# CXX options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-g -O3 -std=c++17 -MD -march=native")

option(NOAVX512 "Disable AVX-512" OFF)
option(NOAVX2 "Disable AVX2" OFF)
if (NOAVX512)
    message("Disabling AVX-512")
    add_compile_options(-mno-avx512f)
endif()
if (NOAVX2)
    message("Disabling AVX2")
    add_compile_options(-mno-avx2)
endif()
if (NOCRT)
    message("Disabling CRT")
    add_compile_options(-DNOCRT=1)
endif()

# main target
set(SRC_FILES src/spiral.cpp src/core.cpp src/constants.cpp src/poly.cpp src/util.cpp src/client.cpp src/testing.cpp)
set(HEADER_FILES include/spiral.h include/core.h include/constants.h include/poly.h include/util.h include/client.h include/testing.h)
add_executable(spiral ${SRC_FILES} ${HEADER_FILES})
target_include_directories(spiral PUBLIC include)
target_include_directories(spiral PUBLIC /home/$ENV{USER}/vcpkg/installed/)

# Must specify:
# -DCMAKE_TOOLCHAIN_FILE=/home/ubuntu/vcpkg/scripts/buildsystems/vcpkg.cmake

find_package(HEXL 1.2.1 CONFIG REQUIRED)
target_link_libraries(spiral PRIVATE HEXL::hexl)

add_compile_definitions(
    TEXP=$(TEXP)
    TEXPRIGHT=$(TEXPRIGHT)
    TCONV=$(TCONV)
    TGSW=$(TGSW)
    QPBITS=$(QPBITS)
    PVALUE=$(PVALUE)
    QNUMFIRST=$(QNUMFIRST)
    QNUMREST=$(QNUMREST)
    OUTN=$(OUTN)
)