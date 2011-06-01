
#include <stdlib.h>
#include <stdio.h>
#include <limits>

#include "data.hpp"
#include "model.hpp"

using namespace std;

struct max_min {
    max_min()
        :max(-numeric_limits<double>::infinity()),min(numeric_limits<double>::infinity()){
    }
    void putVal(double v) {
        max = v > max ? v : max;
        min = v < min ? v : min;
    }
    
    double max;
    double min;
};


/*********************************************************************************
 *  Load Data
 *********************************************************************************/

void Nuclei::updateMinMax() {
    max_min mmx;
    max_min mmy;
    max_min mmz;

	for(Nuclei::iterator i = begin(); i!= end(); ++i) {
        mmx.putVal(i->x);
        mmy.putVal(i->y);
        mmz.putVal(i->z);
	}

    xmax = mmx.max;
    xmin = mmx.min;

    ymax = mmy.max;
    ymin = mmy.min;

    zmax = mmz.max;
    zmin = mmz.min;


}

void Vals::updateMinMax() {
    max_min mmv;

	for(Vals::iterator i = begin(); i!= end(); ++i) {
        mmv.putVal(*i);
	}


    max = mmv.max;
    min = mmv.min;
}



pair<Nuclei,Vals> loadData(const char* filename) {

    Nuclei nuclei;
    Vals vals;
	double x,y,z,v;

	FILE* fp = fopen(filename,"r");

    if(!fp) {
        printf("Could not open %s\n",filename);
        exit(1);
    }


	while(!feof(fp)) {
		if(!fscanf(fp,"%lf %lg %lg %lg \n",&v,&x,&y,&z)) {
            printf("Could not parse input");
            exit(1);
        }
		DataPoint p;
		p.x = x;
		p.y = y;
		p.z = z;


        vals.push_back(v);
        nuclei.push_back(p);
	}
	fclose(fp);

	nuclei.updateMinMax();
	vals.updateMinMax();

    return pair<Nuclei,Vals>(nuclei,vals);
}

double rfloat() {
	return (double)rand()/(double)RAND_MAX;;
}

//Generates fake data by placing spins randomly in the [0,1]^3 cube
//and evaulating a random model. Good for testing how vunerable we are
//to local minima

pair<Nuclei,Vals> fakeData(GaussModel* gaussModel,unsigned long count) {
    Nuclei nuclei;
    Vals vals;

	srand(12345);

	for(unsigned long i=0;i<count;i++) {
		DataPoint p;
		p.x = rfloat();
		p.y = rfloat();
		p.z = rfloat();

		nuclei.push_back(p);
        vals.push_back(gaussModel->eval(p.x,p.y,p.z,1e-3));
	}

	nuclei.updateMinMax();
	vals.updateMinMax();

    return pair<Nuclei,Vals>(nuclei,vals);
}


