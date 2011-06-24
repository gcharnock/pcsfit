
VTKLIBS = -lvtkalglib -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkFiltering -lvtkftgl -lvtkGenericFiltering -lvtkGeovis -lvtkHybrid -lvtkImaging -lvtkInfovis -lvtkIO -lvtkmetaio -lvtkNetCDF -lvtkParallel -lvtkproj4 -lvtkRendering -lvtksqlite -lvtksys -lvtkverdict -lvtkViews -lvtkVolumeRendering -lvtkWidgets

LIBS = -lgsl -lgslcblas -L../Cuba-2.1 -lcuba -lboost_thread

OBJS = main.cpp data.cpp model.cpp #vis.cpp

all:${OBJS}
	g++ -Wall -O3 -ffast-math ${OBJS} -o run  -I../Cuba-2.1 -I/usr/include/vtk-5.4/ ${VTKLIBS}  ${LIBS}  -Wno-deprecated
	./run

