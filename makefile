

all:
	g++ -O3 -ffast-math -o run main.cpp -I../Cuba-2.1 -L../Cuba-2.1 -lcuba -lm -lrt -lgsl -lgslcblas
	./run