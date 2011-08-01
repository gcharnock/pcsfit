
#include "data.hpp"
#include <cmath>
#include <vector>
#include <cassert>
#include <iosfwd>

class Nuclei;
class Vals;

//Okay, we've made the somewhat odd choice of using the curiously
//recuring template patten rather than polymorphisum here which means
//that writing code that is general over the type of model would be
//harder and would require templating. Effectively we are generating
//two or more code paths at compile time rather than one.

//Firstly CRTP has the potential to be slightly more efficent, but
//that probably isn't worth the extra hassel. The real problem is that
//dynamic polymorphisum doesn't seem to gain us much. The only
//function that makes sense to make virtual is eval. At what point can
//we hide the exact model type safely? We have to wait until it gets
//to the minimiser

template<typename Derived>
class ModelBase {
public:
    ModelBase() {
		setBasicParams(1,0,
					   0,0,0,
					   0,0,0);
    }

	void setBasicParams(double _ax,double _rh,
						double x,double y, double z,
						double alpha,double beta,double gamma) {
        ax = _ax;
        rh = _rh; 

        metal.x = x;
        metal.y = y;
        metal.z = z;

        setEulerAngles(alpha,beta,gamma);
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

        for(unsigned long i = 0;i<9;i++) {
            assert(std::isfinite(mat[i]));
        }
    }
	//Evaluate at the points in Nuclei and dump the values in
	//vals. Assume that vals is already allocated to the coorect size.
	void bulkEval(const Nuclei& nuclei,Vals& vals) const {
        for(unsigned long i=0;i<nuclei.size();i++) {
            vals[i] = (static_cast<Derived*>(this))->
				eval(nuclei[i].x,nuclei[i].y,nuclei[i].z);
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

    static std::vector<double> getSmallSteps();
	static std::vector<double> pack(const PointModel& m);
	static PointModel unpack(const std::vector<double>& v);
	static PointModel randomModel(long seed);
};



std::ostream& operator <<(std::ostream& out,const PointModel& m);

class GaussModel : public ModelBase<GaussModel>  {
public:
	GaussModel();
	~GaussModel();

	//This method is thread safe due to being const
	double eval(double x,double y,double z) const;
	double stddev;

    static std::vector<double> getSmallSteps();
	static std::vector<double> pack(const GaussModel& m);
	static GaussModel unpack(const std::vector<double>& v);
	static GaussModel randomModel(long seed);
};



std::ostream& operator <<(std::ostream& out,const GaussModel& m);
