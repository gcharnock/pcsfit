
VTKLIBS = -lvtkalglib -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkFiltering -lvtkftgl -lvtkGenericFiltering -lvtkGeovis -lvtkHybrid -lvtkImaging -lvtkInfovis -lvtkIO -lvtkmetaio -lvtkNetCDF -lvtkParallel -lvtkproj4 -lvtkRendering -lvtksqlite -lvtksys -lvtkverdict -lvtkViews -lvtkVolumeRendering -lvtkWidgets


all:
	g++ -O3 -ffast-math -o run main.cpp data.cpp -I../Cuba-2.1 -L../Cuba-2.1 -lcuba -lm -lrt -lgsl -lgslcblas
	./run

vis: vis.cpp data.cpp
	g++ -g -o vis vis.cpp data.cpp -I/usr/include/vtk-5.4/ ${VTKLIBS} -Wno-deprecated
	./vis