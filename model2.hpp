
#include "data.hpp"

struct PointModel {
    double chi_1,chi_2,chi_xy,chi_xz,chi_yz;
    Vector3 metal;
};

void eval_point(const PointModel* pm,double* value, PointModel* gradient);


