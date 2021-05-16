cmake_minimum_required(VERSION 2.8.9)

project(cuts)

### Policy settings
cmake_policy(SET CMP0054 NEW)   # don't implicitly dereference inside if()


### Process settings
option(BUILD_SHARED_LIBS "Build the shared library" FALSE)
if(BUILD_SHARED_LIBS)
    message("-- Building SHARED libraries")
else()
    message("-- Building STATIC libraries")
endif()


### Resolve optional dependencies

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake") # look for stuff in the /cmake directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/deps/eigen/cmake") # look for stuff in the Eigen/cmake directory
include(UpdateCacheVariable)

# Work with non-standard homebrew installations
# (from ceres build system)
if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
  find_program(HOMEBREW_EXECUTABLE brew)
  mark_as_advanced(FORCE HOMEBREW_EXECUTABLE)
  if (HOMEBREW_EXECUTABLE)
    # Detected a Homebrew install, query for its install prefix.
    execute_process(COMMAND ${HOMEBREW_EXECUTABLE} --prefix
      OUTPUT_VARIABLE HOMEBREW_INSTALL_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Detected Homebrew with install prefix: "
      "${HOMEBREW_INSTALL_PREFIX}, adding to CMake search paths.")
    list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_INSTALL_PREFIX}")
  endif()
endif()


# Look for suitesparse (adapted from ceres build system)
option(SUITESPARSE "Enable SuiteSparse." ON)
if (SUITESPARSE)
  # By default, if SuiteSparse and all dependencies are found, Ceres is
  # built with SuiteSparse support.

  # Check for SuiteSparse and dependencies.
  find_package(SuiteSparse)
  if (SUITESPARSE_FOUND)
    # On Ubuntu the system install of SuiteSparse (v3.4.0) up to at least
    # Ubuntu 13.10 cannot be used to link shared libraries.
    if (BUILD_SHARED_LIBS AND
        SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION)
      message(FATAL_ERROR "You are attempting to build Ceres as a shared "
        "library on Ubuntu using a system package install of SuiteSparse "
        "3.4.0. This package is broken and does not support the "
        "construction of shared libraries (you can still build Ceres as "
        "a static library).  If you wish to build a shared version of Ceres "
        "you should uninstall the system install of SuiteSparse "
        "(libsuitesparse-dev) and perform a source install of SuiteSparse "
        "(we recommend that you use the latest version), "
        "see http://ceres-solver.org/building.html for more information.")
    endif (BUILD_SHARED_LIBS AND
      SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION)

    # By default, if all of SuiteSparse's dependencies are found, Ceres is
    # built with SuiteSparse support.
    message("-- Found SuiteSparse ${SUITESPARSE_VERSION}, "
            "building with SuiteSparse.")
  else (SUITESPARSE_FOUND)
    # Disable use of SuiteSparse if it cannot be found and continue.
    message("-- Did not find all SuiteSparse dependencies, disabling "
      "SuiteSparse support.")
    update_cache_variable(SUITESPARSE OFF)
  endif (SUITESPARSE_FOUND)
else (SUITESPARSE)
  message("-- Building without SuiteSparse.")
endif (SUITESPARSE)

# HACK BY NICK
# The SuiteSparse logic above doesn't look for UMFpack, but we need it.
# This code attempts to find if by assuming that it will be in the same place
# as cholmod
if(SUITESPARSE_FOUND)
    string(REGEX REPLACE "cholmod" "umfpack" UMFPACK_LIBRARY ${CHOLMOD_LIBRARY})
    message("-- Guesstimated umfpack location as: ${UMFPACK_LIBRARY}")
    if(EXISTS ${UMFPACK_LIBRARY})
        list(APPEND SUITESPARSE_LIBRARIES ${UMFPACK_LIBRARY}) 
    else()
        message(WARNING "UMFPack guess failed, so we don't actually have SUITESPARSE support.")
        set(SUITESPARSE_FOUND FALSE)
    endif()
endif()

if(SUITESPARSE_FOUND)
    SET(HAVE_SUITESPARSE TRUE)
    add_definitions(-DHAVE_SUITESPARSE)
    include_directories(${SUITESPARSE_INCLUDE_DIRS})
endif()

### Mosek options
set(USE_MOSEK OFF CACHE BOOL "Use Mosek")
set(MOSEK_DIR "" CACHE STRING "Mosek directory path")
if(USE_MOSEK)
    message("-- Compiling with Mosek")
    SET(HAVE_MOSEK TRUE)
    add_definitions(-DHAVE_MOSEK)
else()
    message("-- Compiling without Mosek")
    SET(HAVE_MOSEK FALSE)
endif()

### Handle windows-specific fixes
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  add_definitions (-DNOMINMAX)  # don't use weird windows built-in min/max
  add_definitions (-D_USE_MATH_DEFINES)  # match unix behavior of constants in cmath
endif()


### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

### Do anything needed for dependencies and bring their stuff in to scope
add_subdirectory(deps)

### Includes from dependencies
include_directories(deps/eigen)
include_directories(deps/glad/include)
include_directories(deps/nanogui/include)
include_directories(deps/nanogui/ext/glfw/include)
include_directories(deps/glm)
include_directories(deps/RectangleBinPack/include)
include_directories(deps/nanort)
include_directories(deps/LBFGS/include)
include_directories(deps/nanoflann)

### Includes from this project
include_directories(core/include)
include_directories(gui/include)

# Includes from the subprojects
include_directories(projects/flatten/core/include)
include_directories(projects/flatten/gui/include)
include_directories(projects/cuts/core/include)
include_directories(projects/cuts/gui/include)

### Compiler options

set( CMAKE_EXPORT_COMPILE_COMMANDS 1 ) # Emit a compile flags file to support completion engines (YouCompleteMe in my case)

### Compiler options
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra -g3")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-private-field -Wno-unused-function -Wno-deprecated-declarations -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-comment")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
  else()
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5 -fopenmp")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRACE_INCLUDES}") # uncomment if you need to track down where something is getting included from
  SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 ")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions (-DNOMINMAX)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()


### Recurse in to each of the subprojects
add_subdirectory(core)
add_subdirectory(projects)
add_subdirectory(gui)
#add_subdirectory(plugins)
