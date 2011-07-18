
#include "data.hpp"
#include <cmath>
#include <vector>

class Nuclei;
class Vals;

template<typename Derived>
class ModelBase {
public:
    ModelBase() {
        ax = 1;
        rh = 0; 

        metal.x = 0;
        metal.y = 0;
        metal.z = 0;

        setEulerAngles(0,0,0);
    }

    
    void setEulerAngles(double _angle_x,double _angle_y,double _angle_z) {
        angle_x=_angle_x;
        angle_y=_angle_y;
        angle_z=_angle_z;

        //From the matrix and quaternion FAQ
        double A       = std::cos(angle_x);
        double B       = std::sin(angle_x);
        double C       = std::cos(angle_y);
        double D       = std::sin(angle_y);
        double E       = std::cos(angle_z);
        double F       = std::sin(angle_z);
        double AD      = A * D;
        double BD      = B * D;
        mat[0]  =   C * E;
        mat[1]  =  -C * F;
        mat[2]  =   D;

        mat[3]  =  BD * E + A * F;
        mat[4]  = -BD * F + A * E;
        mat[5]  =  -B * C;
        mat[6]  = -AD * E + B * F;
        mat[7]  =  AD * F + B * E;
        mat[8] =   A * C;
    }
	//Evaluate at the points in Nuclei and dump the values in
	//vals. Assume that vals is already allocated to the coorect size.
	void bulkEval(const Nuclei& nuclei,Vals& vals) const {
        for(unsigned long i=0;i<nuclei.size();i++) {
            vals[i] = (static_cast<Derived*>(this))->eval(nuclei[i].x,nuclei[i].y,nuclei[i].z);
        }
        vals.updateMinMax();
    }

    Vector3 metal;
    
    double ax;
    double rh;

    double angle_x,angle_y,angle_z;

    double mat[9];

};

class PointModel : public ModelBase<PointModel> {
public:
    PointModel();
    ~PointModel();

    double eval(double x, double y, double z) const;
    
};
std::vector<double> packPointModel(const PointModel& m);
PointModel unpackPointModel(const std::vector<double>& v);

class GaussModel : public ModelBase<GaussModel>  {
public:
	GaussModel();
	~GaussModel();

	//This method is thread safe due to being const
	double eval(double x,double y,double z) const;

	double exponant;

	double cube_x_min, cube_x_max;
	double cube_y_min, cube_y_max;
	double cube_z_min, cube_z_max;

	//3x3 rotation matrix
	double mat[9];
};
std::vector<double> packGaussModel(const GaussModel& m);
GaussModel unpackGaussModel(const std::vector<double>& v);
