
#include <vector>
#include <utility>

struct DataPoint {
	double x;
	double y;
	double z;
};

struct Nuclei : std::vector<DataPoint> {
    double xmax;
    double xmin;
    
    double ymax;
    double ymin;

    double zmax;
    double zmin;
};

struct Vals : std::vector<double> {    
    double max;
    double min;
};

typedef std::pair<Nuclei,Vals> pairNucVals;

pairNucVals loadData(const char* filename);
pairNucVals fakeData(unsigned long count);
