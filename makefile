
VTKLIBS = -lvtkalglib -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkFiltering -lvtkftgl -lvtkGenericFiltering -lvtkGeovis -lvtkHybrid -lvtkImaging -lvtkInfovis -lvtkIO -lvtkmetaio -lvtkNetCDF -lvtkParallel -lvtkproj4 -lvtkRendering -lvtksqlite -lvtksys -lvtkverdict -lvtkViews -lvtkVolumeRendering -lvtkWidgets

LIBS = -lgsl -lgslcblas -L../Cuba-2.1 -lcuba -L. -lboost_thread-gcc44-mt-d-1_46_1 -lpthread

OBJS = main.cpp data.cpp model.cpp vis.cpp

INCS = -I../Cuba-2.1 -I/usr/include/vtk-5.4/ -I../../src/boost_1_46_1

#CFLAGS = -O3 -ffast-math
CFLAGS = -g

all:${OBJS}
	g++ -Wall ${CFLAGS} ${OBJS} -o run ${INCS}  ${VTKLIBS}  ${LIBS}  -Wno-deprecated
	./run

threads: threads.cpp
	g++ -Wall threads.cpp -o threads ${LIBS}
