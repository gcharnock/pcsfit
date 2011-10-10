
cmake_minimum_required(VERSION 2.8)
project(pcsfit)

#set(CMAKE_VERBOSE_MAKEFILE TRUE)

#Debug build, can use with gdb
#set(CMAKE_BUILD_TYPE Debug)

#Compile a release build
set(CMAKE_BUILD_TYPE Release)


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
  -lboost_program_options_link
  -lpthread)


if(IS_DIRECTORY /usr/include/vtk-5.4)
  set(vtkDir /usr/include/vtk-5.4/)
else()
  set(vtkDir /usr/include/vtk-5.2/)
endif()


include_directories(../Cuba-2.1 ${vtkDir}  $ENV{HOME}/src/boost_1_46_1)


set(CMAKE_CXX_FLAGS         "-Wno-deprecated -Wall -fno-strict-aliasing")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -ggdb")
set(cmake_cxx_flags_release "-O3 -ffast-math")
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS})


set(sources  model.cpp model2.cpp fit.cpp  data.cpp  vis.cpp tests.cpp)


add_executable(pcsfit main.cpp  ${sources})
target_link_libraries(pcsfit ${LIBS} ${VTKLIBS})

