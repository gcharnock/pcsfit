
class Nuclei;
class Vals;

#include "data.hpp"

class GaussModel {
public:
	GaussModel();

	//This method is thread safe due to being const
	double eval(double x,double y,double z,double epsAbs) const;

	//Evaluate at the points in Nuclei and dump the values in
	//vals. Assume that vals is already allocated to the coorect size.
	void bulkEval(const Nuclei& nuclei,Vals& vals) const;

	Vector3 metal;

	double ax;
	double rh;

	//Set the euler angles
	void setEulerAngles(double angle_x,double angle_y,double angle_z);

	double exponant;

	double cube_x_min, cube_x_max;
	double cube_y_min, cube_y_max;
	double cube_z_min, cube_z_max;

	//Todo, needs to be a function of exponant
	double cutoff2;

	//3x3 rotation matrix
	double mat[9];
};
