cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()


# Collect all .cpp files in the subdirectories
file(GLOB_RECURSE FLATTEN_SRCS *.cpp)
file(GLOB_RECURSE FLATTEN_HEADS *.h *.hpp *.ipp)

# Create a single library for the project
add_library(
        flatten-core
        ${FLATTEN_SRCS}
	${FLATTEN_HEADS}
        )


target_link_libraries(flatten-core core)
if (HAVE_MOSEK)
	list(APPEND MOSEK_LIBS ${MOSEK_DIR}/libiomp5.dylib)
	list(APPEND MOSEK_LIBS ${MOSEK_DIR}/libmosek64.7.1.dylib)
	list(APPEND MOSEK_LIBS ${MOSEK_DIR}/libmosek64.dylib)
	target_link_libraries(flatten-core ${MOSEK_LIBS})
	file(COPY ${MOSEK_LIBS} DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Release)
	file(COPY ${MOSEK_LIBS} DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Debug)
endif()
