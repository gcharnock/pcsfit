

#include "maths.hpp"
#include "foreach.hpp"

#include <boost/array.hpp>

using namespace std;
using namespace boost;

double magic_polynomial(double t) {
    double t2 =  t*t;
    double t3 =  t2*t;
    double t4 =  t3*t;
    double t5 =  t4*t;
    double t6 =  t5*t;
    double t7 =  t6*t;
    double t8 =  t7*t;
    double t9 =  t8*t;
    double t10 = t9*t;
    double t11 = t10*t;

    return 3879876*t11*
        (+      t10 /21.0
         -      t9  /2.0
         +45  * t8  /19.0
         -20  * t7  /3.0
         +210 * t6  /17.0
         -63  * t5  /4.0
         +14  * t4
         -60  * t3  /7.0
         +45  * t2  /13.0
         -5   *t    /6.0
         +1         /11.0);
}


double bump(double t)  {
    return abs(t) < 1 ? exp(-1/(1-t*t)) : 0;
}

//Since we need to pass this function to a C api it can't actually be
//a member function. We'll pass the this pointer explicitly to fake a
//member
//Is [c,d] a subset of [a,b]?
bool subset(double a,double b,double c,double d) {
    assert(a<b);
    assert(c<d);
    return (c<a && c<b) || (d<a && d<b);
}
bool intersect(double a,double b,double c,double d) {
    return subset(a,b,c,d) || subset(c,d,a,b);
}


Matrix3d MakeMatrix3d(double a00, double a01, double a02,
                      double a10, double a11, double a12,
                      double a20, double a21, double a22) {
    Matrix3d m;
    m(0,0)=a00;             m(0,1)=a01;             m(0,2)=a02;
    m(1,0)=a10;             m(1,1)=a11;             m(1,2)=a12;
    m(2,0)=a20;             m(2,1)=a21;             m(2,2)=a22;

    return m;
}

