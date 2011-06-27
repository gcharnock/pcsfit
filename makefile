
VTKLIBS = -lvtkalglib -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkFiltering -lvtkftgl -lvtkGenericFiltering -lvtkGeovis -lvtkHybrid -lvtkImaging -lvtkInfovis -lvtkIO -lvtkmetaio -lvtkNetCDF -lvtkParallel -lvtkproj4 -lvtkRendering -lvtksqlite -lvtksys -lvtkverdict -lvtkViews -lvtkVolumeRendering -lvtkWidgets

LIBS = -lgsl -lgslcblas -L../Cuba-2.1 -lcuba -lboost_thread-mt
INCS = -I../Cuba-2.1
CFLAGS = -g #-O3 -ffast-math 

all:
	g++ ${CFLAGS} -o run main.cpp data.cpp model.cpp ${INCS} ${LIBS} 
	./run

vis: vis.cpp data.cpp model.cpp threads.cpp
	g++ -g -o vis vis.cpp data.cpp model.cpp threads.cpp -I../Cuba-2.1 -I/usr/include/vtk-5.4/ ${VTKLIBS}  ${LIBS} -Wno-deprecated
