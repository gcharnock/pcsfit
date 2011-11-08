
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


void loadData(const std::string& filename,Dataset* dataset) {
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
        if(v == 0.0) {
            continue;
        }
        dataset->vals.push_back(v);
        dataset->nuclei.push_back(Vector3(x,y,z));
	}
	fclose(fp);

	dataset->nuclei.updateMinMax();
	dataset->vals.updateMinMax();
}

void toSketch3d(ostream& out,Dataset* dataset) {
	out.precision(5);
	out << fixed;


    if (dataset->nuclei.size() == 0) {
        cout << "WARNING: dataset contains no spins" << endl;
    }

	dataset->nuclei.updateMinMax();
	dataset->vals.updateMinMax();

    out << "def nsegs 10"  << endl;

    double xmin = dataset->nuclei.xmin / 3;
    double xmax = dataset->nuclei.xmax / 3;

    double ymin = dataset->nuclei.ymin / 3;
    double ymax = dataset->nuclei.ymax / 3;

    double zmin = dataset->nuclei.zmin / 3;
    double zmax = dataset->nuclei.zmax / 3;

    out << "def scene {" << endl;

    for(unsigned long i = 0; i < dataset->nuclei.size(); i++) {
        double x = dataset->nuclei[i].x / 3;
        double y = dataset->nuclei[i].y / 3;
        double z = dataset->nuclei[i].z / 3;

        double absMax = abs(dataset->vals.max) > abs(dataset->vals.min) ? abs(dataset->vals.max) : abs(dataset->vals.min);

        double v = dataset->vals[i] / absMax; // V has a range of [-1,1]
        double r = pow(abs(v)+0.0001,1.0/3);

        double intensity = pow(abs(v),1.0/3)*100;
        string redblue   = (v >= 0 ? "red" : "blue");

        out << "put{translate([" << x << "," << y << "," << z << "])}{" << endl;
        out << "    sweep[draw=black!20,fill=" << redblue << "!" << intensity << "!black,style=very thin]" << endl;
        out << "        {nsegs, rotate(360 / nsegs, [0,1,0]) }" << endl;
        out << "         sweep {nsegs, rotate(360 / nsegs, [1,0,0]) }" << endl;
        out << "               (0,0," << r << ")" << endl;
        out << "}" << endl;
    }
    
    //Coordinate axies

    //X axis
    out << "def x (" << xmax << "," << ymin << "," << zmin << ")" << endl;
    out << "line[draw=black,arrows=->] (" << xmin << "," << ymin << "," << zmin << ") (x)" << endl;
    //Y axis
    out << "def y (" << xmin << "," << ymax << "," << zmin << ")" << endl;
    out << "line[draw=black,arrows=->] (" << xmin << "," << ymin << "," << zmin << ") (y)" << endl;
    //Z axis
    out << "def z (" << xmin << "," << ymin << "," << zmax << ")" << endl;
    out << "line[draw=black,arrows=->] (" << xmin << "," << ymin << "," << zmin << ") (z)" << endl;

    out << "special |\\node at #1 {$x$};\\node at #2 {$y$};\\node at #3 {$z$};|(x)(y)(z)" << endl;


    out << "}" << endl;

    out << "put { view((" << xmax << "," << ymax << "," << zmax << "),"
        << "(" << (xmax+xmin)/2 << "," << (ymax+ymin)/2 << "," << (zmax+zmin)/2 << ")"
        << ") } {scene}" << endl;
    out << "global{language tikz}" << endl;
}

