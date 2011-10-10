
#ifndef __DATA__H__
#define __DATA__H__

#include <vector>
#include <string>
#include <utility>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

typedef boost::mt19937 PRNG;
typedef boost::uniform_01<> RandomDist;

struct GaussModel;

struct Vector3 {
	Vector3(double _x,double _y,double _z)
		: x(_x),y(_y),z(_z){}
	Vector3() {}
	double x;
	double y;
	double z;
};

struct Nuclei : std::vector<Vector3> {
	void updateMinMax();
    double xmax;
    double xmin;
    
    double ymax;
    double ymin;

    double zmax;
    double zmin;
};

struct Vals : std::vector<double> {    
	void updateMinMax();

    double max;
    double min;
};

struct Dataset {
    Nuclei nuclei;
    Vals vals;
};

typedef std::pair<Nuclei,Vals> pairNucVals;

void loadData(const std::string& filename,Dataset* dataset);



#endif
