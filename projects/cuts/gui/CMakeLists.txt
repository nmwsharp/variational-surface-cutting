cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()


# Collect all .cpp files in the subdirectories
file(GLOB_RECURSE CUTS_GUI_SRCS *.cpp)

# Create a single library for the project
add_library(
        cuts-gui
        ${CUTS_GUI_SRCS}
        )

target_link_libraries(cuts-gui viewer cuts-core core)
