
VTKLIBS = -lvtkalglib -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkFiltering -lvtkftgl -lvtkGenericFiltering -lvtkGeovis -lvtkHybrid -lvtkImaging -lvtkInfovis -lvtkIO -lvtkmetaio -lvtkNetCDF -lvtkParallel -lvtkproj4 -lvtkRendering -lvtksqlite -lvtksys -lvtkverdict -lvtkViews -lvtkVolumeRendering -lvtkWidgets

LIBS = -lgsl -lgslcblas -L../Cuba-2.1 -lcuba -L. -lboost_thread_link -lpthread

OBJS = main.o data.o model.o vis.o

CXX = g++

HEADERS = minimiser.hpp

INCS = -I../Cuba-2.1 -I/usr/include/vtk-5.4/ -I../../src/boost_1_46_1

#CFLAGS = -O3 -ffast-math -Wno-deprecated
CFLAGS = -g -Wno-deprecated

all:${OBJS}
	${CXX} -o run ${OBJS} ${VTKLIBS}  ${LIBS}  

.cpp.o:
	${CXX} -c -Wall $(INCS) $(CFLAGS) $<


clean:
	rm ${OBJS}