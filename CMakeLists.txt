cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(trish VERSION 0.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
endif()


find_package(CGAL 5.2 REQUIRED)


add_executable(trish a.cc b.cc trish.cc)
target_link_libraries(trish CGAL)

