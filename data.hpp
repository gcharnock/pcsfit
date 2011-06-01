
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

std::pair<Nuclei,Vals> loadData(const char* filename);
std::pair<Nuclei,Vals> fakeData(unsigned long count);
