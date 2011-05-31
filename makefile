
VTKLIBS = -lvtkalglib -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkFiltering -lvtkftgl -lvtkGenericFiltering -lvtkGeovis -lvtkHybrid -lvtkImaging -lvtkInfovis -lvtkIO -lvtkmetaio -lvtkNetCDF -lvtkParallel -lvtkproj4 -lvtkQtChart -lvtkRendering -lvtksqlite -lvtksys -lvtkverdict -lvtkViews -lvtkVolumeRendering -lvtkWidgets


all:
	g++ -O3 -ffast-math -o run main.cpp -I../Cuba-2.1 -L../Cuba-2.1 -lcuba -lm -lrt -lgsl -lgslcblas
	./run

vis: vis.cpp
	g++ -g -o vis vis.cpp -I/usr/include/vtk-5.4/ ${VTKLIBS}
	./vis