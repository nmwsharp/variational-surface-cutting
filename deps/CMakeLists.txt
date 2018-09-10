cmake_minimum_required(VERSION 2.8.9)

### Includes
include_directories(eigen)
include_directories(glad/include)
include_directories(glm)

### Recurse to subdirecties with source
add_subdirectory(glad)


### Nanogui

# Don't build nanogui examples or python bindings (python takes forever)
set(NANOGUI_BUILD_PYTHON OFF CACHE FALSE "don't build python")
set(NANOGUI_BUILD_EXAMPLE OFF CACHE FALSE "don't build examples")

# Nanogui's cmake script sets BUILD_SHARED_LIBS to the value of NANOGUI_BUILD_SHARED, we want the inverse behavior
# (note that this effectively makes the NANOGUI_BUILD_SHARED variable useless
set(NANOGUI_BUILD_SHARED ${BUILD_SHARED_LIBS} CACHE BOOL " " FORCE)

add_subdirectory(nanogui)

# Pass nanogui variables up to parent 
set(NANOGUI_EXTRA_LIBS ${NANOGUI_EXTRA_LIBS} PARENT_SCOPE)
set(NANOGUI_EXTRA_DEFS ${NANOGUI_EXTRA_DEFS} PARENT_SCOPE)
set(NANOGUI_EXTRA_INCS ${NANOGUI_EXTRA_INCS} PARENT_SCOPE)


# Rectangle bin pack
add_subdirectory(RectangleBinPack)