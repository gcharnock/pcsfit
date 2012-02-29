
cmake_minimum_required(VERSION 2.8)
project(pcsfit)
find_package(PkgConfig)

#set(CMAKE_VERBOSE_MAKEFILE TRUE)

#Debug build, can use with gdb
set(CMAKE_BUILD_TYPE Debug)

#Compile a release build
#set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_CXX_COMPILER "g++-4.5")
#set(CMAKE_CXX_LINKER "g++-4.5")

#set(VTKLIBS
#  vtkCommon
#  vtkftgl
#  vtkRendering)

set(LIBS  
  gsl
  gslcblas
  -L../Cuba-2.1
  cuba
  -L.
  -L${CMAKE_CURRENT_BINARY_DIR}
  -lboost_thread_link
  -lboost_program_options_link
  -lpthread
)

pkg_search_module(SIGC REQUIRED sigc++-2.0)
include_directories(SYSTEM ${SIGC_INCLUDE_DIRS})


#if(IS_DIRECTORY /usr/include/vtk-5.4)
#  set(vtkDir /usr/include/vtk-5.4/)
#else()
#  set(vtkDir /usr/include/vtk-5.2/)
#endif()


include_directories(
  ../Cuba-2.1 ${vtkDir}
  $ENV{HOME}/src/boost
  ../libspinxml/
)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR})


set(CMAKE_CXX_FLAGS         "-Wno-deprecated -Wall -fno-strict-aliasing")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -ggdb")
set(cmake_cxx_flags_release "-O3 -ffast-math")
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS})

set(spinxmlSources 
  ../libspinxml/spinxml/orientation.cpp
  ../libspinxml/spinxml/interaction.cpp
  ../libspinxml/spinxml/panic.cpp
  ../libspinxml/spinxml/basic_math.cpp
  ../libspinxml/spinxml/unit.cpp
  ../libspinxml/spinxml/frame.cpp
  )
set(sources model2.cpp fit.cpp  data.cpp  tests.cpp pointdev.cpp maths.cpp ${spinxmlSources})


add_executable(pcsfit main.cpp  ${sources})
target_link_libraries(pcsfit ${LIBS} ${VTKLIBS} sigc-2.0)


