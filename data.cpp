
#include <stdlib.h>
#include <stdio.h>
#include <limits>

#include "data.hpp"

using namespace std;

/*********************************************************************************
 *  Load Data
 *********************************************************************************/

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

pair<Nuclei,Vals> loadData(const char* filename) {
    Nuclei nuclei;
    Vals vals;
	double x,y,z,v;

	FILE* fp = fopen(filename,"r");

    if(!fp) {
        printf("Could not open %s\n",filename);
        exit(1);
    }

    max_min mmx;
    max_min mmy;
    max_min mmz;
    max_min mmv;

	while(!feof(fp)) {
		if(!fscanf(fp,"%lf %lg %lg %lg \n",&v,&x,&y,&z)) {
            printf("Could not parse input");
            exit(1);
        }
		DataPoint p;
		p.x = x;
		p.y = y;
		p.z = z;

        mmx.putVal(x);
        mmy.putVal(y);
        mmz.putVal(z);
        mmv.putVal(v);

        vals.push_back(v);
        nuclei.push_back(p);
	}
	fclose(fp);

    nuclei.xmax = mmx.max;
    nuclei.xmin = mmx.min;

    nuclei.ymax = mmy.max;
    nuclei.ymin = mmy.min;

    nuclei.zmax = mmz.max;
    nuclei.zmin = mmz.min;

    vals.max = mmv.max;
    vals.min = mmv.min;

    return pair<Nuclei,Vals>(nuclei,vals);
}

double rfloat() {
	return (double)rand()/(double)RAND_MAX;;
}

//Generates fake data by placing spins randomly in the [0,1]^3 cube
//and evaulating a random model. Good for testing how vunerable we are
//to local minima

pair<Nuclei,Vals> fakeData(unsigned long count) {
    Nuclei nuclei;
    Vals vals;

	srand(12345);

	for(unsigned long i=0;i<count;i++) {
		DataPoint p;
		p.x = rfloat();
		p.y = rfloat();
		p.z = rfloat();

		nuclei.push_back(p);
        vals.push_back(0); /* TODO */
	}
    return pair<Nuclei,Vals>(nuclei,vals);
}


