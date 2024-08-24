cmake_minimum_required(VERSION 3.30)

set(ENV{CC} /usr/local/bin/gcc)
set(ENV{CXX} /usr/local/bin/g++)
project(CFD
    VERSION 0.1.0
    LANGUAGES CXX
)

# Specify the C standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)


include(FetchContent)
Set(FETCHCONTENT_QUIET FALSE)


# Armadillo (for 2D problems)
if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    message(STATUS "Detected macOS")
    # Enable detection of OpenBLAS on macOS
    option(ALLOW_OPENBLAS_MACOS "Enable detection of OpenBLAS on macOS" ON)
    option(ALLOW_BLAS_LAPACK_MACOS "Enable detection of standard BLAS/LAPACK on macOS" ON)

    if(ALLOW_OPENBLAS_MACOS)
        message(STATUS "ALLOW_OPENBLAS_MACOS is enabled")
        set(ALLOW_OPENBLAS_MACOS ON)
    endif()

    if(ALLOW_BLAS_LAPACK_MACOS)
        message(STATUS "ALLOW_BLAS_LAPACK_MACOS is enabled")
        set(ALLOW_BLAS_LAPACK_MACOS ON)
    endif()
else()
    message(STATUS "Operating system is not macOS. Skipping macOS-specific options.")
endif()

FetchContent_Declare(
    armadillo
    GIT_REPOSITORY  https://gitlab.com/conradsnicta/armadillo-code.git
    GIT_TAG         "14.0.2"
    GIT_CONFIG      "advice.detachedhead=false"
    GIT_SHALLOW     TRUE
    GIT_PROGRESS    TRUE
)
FetchContent_MakeAvailable(armadillo)


add_subdirectory(utils)
add_subdirectory(src)
