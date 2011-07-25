
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



pair<Nuclei,Vals> loadData(const string& filename) {

    Nuclei nuclei;
    Vals vals;
	double x,y,z,v;

	FILE* fp = fopen(filename.c_str(),"r");

    if(!fp) {
        printf("Could not open %s\n",filename.c_str());
        exit(1);
    }


	while(!feof(fp)) {
		if(!fscanf(fp,"%lf %lg %lg %lg \n",&v,&x,&y,&z)) {
            printf("Could not parse input");
            exit(1);
        }
		Vector3 p;
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

double rdouble() {
	return (double)rand()/(double)RAND_MAX;;
}

