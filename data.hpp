
#ifndef __DATA__H__
#define __DATA__H__

#include <vector>
#include <string>
#include <utility>
#include <ostream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

#include "maths.hpp"

typedef boost::mt19937 PRNG;
typedef boost::uniform_01<> RandomDist;

struct GaussModel;

struct Nuclei : std::vector<Vec3d> {
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

void toSketch3d(std::ostream& out,Dataset* dataset);

#endif
