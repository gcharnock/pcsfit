
cmake_minimum_required(VERSION 2.8)
project(pcsfit)

#set(CMAKE_VERBOSE_MAKEFILE TRUE)

#Debug build, can use with gdb
set(CMAKE_BUILD_TYPE Debug)

#Compile a release build
#set(CMAKE_BUILD_TYPE Release)


set(VTKLIBS
  vtkCommon
  vtkftgl
  vtkRendering)

set(LIBS  
  gsl
  gslcblas
  -L../Cuba-2.1
  cuba
  -L.
  -lboost_thread_link
  -lboost_program_options
  -lpthread)

include_directories(../Cuba-2.1 /usr/include/vtk-5.4/ ../../src/boost_1_46_1)


set(CMAKE_CXX_FLAGS         "-std=c++0x -Wno-deprecated -Wall")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -ggdb")
set(cmake_cxx_flags_release "-O3 -ffast-math")


set(sources  model.cpp  data.cpp  vis.cpp)


add_executable(gauss main.cpp  ${sources})
target_link_libraries(gauss ${LIBS} ${VTKLIBS})

add_executable(point point.cpp ${sources})
target_link_libraries(point ${LIBS} ${VTKLIBS})

