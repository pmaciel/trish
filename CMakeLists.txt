cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(trish VERSION 0.0.0 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
endif()


find_package(CGAL 5.2 REQUIRED)
find_package(Qhull REQUIRED)


add_executable(trish a.cc b.cc c.cc trish.cc)
target_link_libraries(trish CGAL qhull_r qhullcpp)
target_link_directories(trish PRIVATE "${Qhull_DIR}/../..")


add_executable(stripack
    stripack/addnod.f90
    stripack/bdyadd.f90
    stripack/covsph.f90
    stripack/insert.f90
    stripack/intadd.f90
    stripack/jrand.f90
    stripack/left.f90
    stripack/lstptr.f90
    stripack/store.f90
    stripack/swap.f90
    stripack/swptst.f90
    stripack/trfind.f90
    stripack/stritest.f90)

