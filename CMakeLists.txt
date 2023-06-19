cmake_minimum_required(VERSION 3.24)
project(Imporve_SeamCarving)

set(CMAKE_CXX_STANDARD 14)

#set(DCMAKE_BUILD_TYPE=Release)

file(GLOB src "${PROJECT_SOURCE_DIR}/*.cpp")

add_executable(SeamCarving ${src})

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})

target_link_libraries(SeamCarving ${OpenCV_LIBS})