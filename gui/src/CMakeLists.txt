cmake_minimum_required(VERSION 2.8.9)

# Collect all .cpp files in the subdirectories
#file(GLOB_RECURSE VIEWER_SRCS *.cpp)

# To keep compile times reasonable here we will create a library with all of the shared
# gui-related resources, as well as a library for each subproject.


### Create a library for shared gui resources
add_library(
	viewer
    camera.cpp colormaps.cpp gl_utils.cpp scene_object.cpp scene_mesh.cpp scene_points.cpp scene_lines.cpp scene_plane.cpp scene_vectors.cpp viewer.cpp project_gui.cpp image.cpp
	)

target_link_libraries(viewer core ${GLFW_LIBRARIES})

if(APPLE)
    target_link_libraries(viewer nanogui ${NANOGUI_EXTRA_LIBS}) # don't need extension loader on OSX
elseif(WIN32)
    target_link_libraries(viewer nanogui ${NANOGUI_EXTRA_LIBS}) # windows uses glad from nanogui
else()
    target_link_libraries(viewer nanogui glad ${NANOGUI_EXTRA_LIBS})
endif()



### Create libraries for individual projects
set(SUB_GUI_LIBRARIES "") # variable used to accumulate target names


LIST(APPEND SUB_GUI_LIBRARIES "cuts-gui")


### Create an executable for the gui
add_executable(cuts gui.cpp main.cpp)
target_link_libraries(cuts viewer ${SUB_GUI_LIBRARIES})
