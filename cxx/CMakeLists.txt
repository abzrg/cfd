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

add_subdirectory(utils)
add_subdirectory(src)
