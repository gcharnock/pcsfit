
struct IntegralBounds;
typedef int (*IntegrandF) (const double xx[],double ff[], int ncomp, void *data);

struct IntegralBounds {
    double xmax;
    double xmin;
               
    double ymax;
    double ymin;
               
    double zmax;
    double zmin;

    void* data;
	IntegrandF integrand;
};


void cuhreIntegrate(IntegrandF f,IntegralBounds* bounds,unsigned long ncomp,double* integral,void* data);
