add_library(utils STATIC
    linspace.cpp
    io.cpp
)
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(utils
    PUBLIC armadillo
)
