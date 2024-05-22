cmake_minimum_required(VERSION 3.15...3.26)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -fopenmp")
find_package(GSL REQUIRED)

file(GLOB SOURCES "src/*.cpp")
add_executable(integrate main.cpp ${SOURCES})
target_include_directories(integrate PRIVATE ${GSL_INCLUDE_DIRS} include)
target_link_libraries(integrate PRIVATE ${GSL_LIBRARIES})
