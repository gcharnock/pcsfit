
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

typedef std::pair<Nuclei,Vals> pairNucVals;

pairNucVals loadData(const std::string& filename);


//Generates fake data by placing spins randomly in the [0,50]^3 cube
//and evaulating a random model. Good for testing how vunerable we are
//to local minima

template<typename Model>
pairNucVals fakeData(long seed,const Model& model,unsigned long count) {
    Nuclei nuclei;
    Vals vals;
	PRNG prng(seed);
	RandomDist rand;

	for(unsigned long i=0;i<count;i++) {
		Vector3 p;
		p.x = 50*rand(prng);
		p.y = 50*rand(prng);
		p.z = 50*rand(prng);

		nuclei.push_back(p);
        vals.push_back(model.eval(p.x,p.y,p.z));
	}

	nuclei.updateMinMax();
	vals.updateMinMax();

    return std::pair<Nuclei,Vals>(nuclei,vals);
}

#endif
