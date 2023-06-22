cmake_minimum_required(VERSION 3.24)
project(Imporve_SeamCarving)

set(CMAKE_CXX_STANDARD 14)

file(GLOB src "${PROJECT_SOURCE_DIR}/*.cpp")

add_executable(SeamCarving ${src} )

# glad
include_directories(./glad/include)

# lsd
#include_directories(./lsd/lsd.h)

set(LOCAL_H /opt/homebrew/include)
include_directories(${LOCAL_H})

# opencv
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

# Eigen
include_directories(/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3)

set(GLEW_H /opt/homebrew/Cellar/glew/2.2.0_1/include/GL)
set(GLFW_H /opt/homebrew/Cellar/glfw/3.3.8/include/GLFW)
link_directories(${GLEW_H} ${GLFW_H})

set(GLEW_LINK /opt/homebrew/Cellar/glew/2.2.0_1/lib/libGLEW.2.2.dylib)
set(GLFW_LINK /opt/homebrew/Cellar/glfw/3.3.8/lib/libglfw.3.dylib)
link_libraries(${OPENGL} ${GLEW_LINK} ${GLFW_LINK})

target_link_libraries(SeamCarving ${OpenCV_LIBS})

if (APPLE)
    target_link_libraries(SeamCarving "-framework OpenGL")
    target_link_libraries(SeamCarving "-framework GLUT")
endif ()