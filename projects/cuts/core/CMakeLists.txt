cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()


# Collect all .cpp files in the subdirectories
file(GLOB_RECURSE CUTS_SRCS *.cpp)

# Create a single library for the project
add_library(
        cuts-core
        ${CUTS_SRCS}
        )

target_link_libraries(cuts-core core rectangle-bin-pack flatten-core)
