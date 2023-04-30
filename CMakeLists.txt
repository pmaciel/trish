cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(trish VERSION 0.0.0 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
endif()


find_package(CGAL 5.2 REQUIRED)


add_executable(trish a.cc b.cc trish.cc)
target_link_libraries(trish CGAL)


add_executable(stripack
    stripack/addnod.f90
    stripack/areas.f90
    stripack/bdyadd.f90
    stripack/bnodes.f90
    stripack/circum.f90
    stripack/covsph.f90
    stripack/crlist.f90
    stripack/delarc.f90
    stripack/delnb.f90
    stripack/delnod.f90
    stripack/edge.f90
    stripack/getnp.f90
    stripack/insert.f90
    stripack/inside.f90
    stripack/intadd.f90
    stripack/intrsc.f90
    stripack/jrand.f90
    stripack/left.f90
    stripack/lstptr.f90
    stripack/nbcnt.f90
    stripack/nearnd.f90
    stripack/optim.f90
    stripack/scoord.f90
    stripack/store.f90
    stripack/swap.f90
    stripack/swptst.f90
    stripack/trans.f90
    stripack/trfind.f90
    stripack/trlist.f90
    stripack/trlprt.f90
    stripack/trmesh.f90
    stripack/trplot.f90
    stripack/trprnt.f90
    stripack/vrplot.f90

    stripack/stritest.f90)

