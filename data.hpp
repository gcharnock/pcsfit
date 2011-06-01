
#include <vector>
#include <utility>

struct GaussModel;

struct DataPoint {
	double x;
	double y;
	double z;
};

struct Nuclei : std::vector<DataPoint> {
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

pairNucVals loadData(const char* filename);
pairNucVals fakeData(GaussModel* gaussModel,unsigned long count);
