
#ifndef MATHS_HPP
#define MATHS_HPP

#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <map>
#include <cmath>
#include <boost/shared_array.hpp>
#include <Eigen/Dense>

//Of course, you shouldn't use namespaces in header files, but we'll
//probably be alright in this case.
using namespace Eigen;

typedef unsigned long ulong;

template<typename T>
struct Vec3 {
	Vec3() {}
	Vec3(T x,T y,T z) {
        mData[0] = x;
        mData[1] = y;
        mData[2] = z;
    }
    Vec3(const T* array) {
        memcpy(mData,array,sizeof(T)*3);
    }

    Vec3(const Vec3& other) {
        memcpy(mData,other.mData,3*sizeof(T));
    }

    Vec3 operator+(const Vec3& other) {
        return Vec3(mData[0] + other.mData[0],
                    mData[1] + other.mData[1],
                    mData[2] + other.mData[2]);
    }

    Vec3 operator-(const Vec3& other) {
        return Vec3(mData[0] - other.mData[0],
                    mData[1] - other.mData[1],
                    mData[2] - other.mData[2]);
    }

    T dot(Vec3 other) {
        return mData[0]*other.mData[0] +
            mData[1]*other.mData[1] +
            mData[2]*other.mData[2];
    }

    T r2() {
        return
            mData[0]*mData[0] +
            mData[1]*mData[1] +
            mData[2]*mData[2];
    }

    Vec3 operator*(T s) {
        return Vec3(mData[0] * s,
                    mData[1] * s,
                    mData[2] * s);
    }


    const Vec3& operator=(const Vec3& other) {
        memcpy(mData,other.mData,3*sizeof(T));
        return *this;
    }

    T operator[](ulong i) {
        assert(i < 3);
        return mData[i];
    }
    void set(ulong i, T v) {
        assert(i < 3);
        mData[i] = v;
    }

    //We need a total ordering so we can put these in a map
    bool operator<(const Vec3& other) const {
        for(ulong i = 0; i < 3; i++) {
            T d1 = mData[i]-other.mData[i];
            if(d1 != 0) {
                return mData[i] < other.mData[i];
            }
        }
        return false;
    }

    const T& x() const {return mData[0];}
    const T& y() const {return mData[1];}
    const T& z() const {return mData[2];}

    T& x() {return mData[0];}
    T& y() {return mData[1];}
    T& z() {return mData[2];}
private:
    T mData[3];
};

typedef Vec3<long> Vec3l;
typedef Vec3<ulong> Vec3ul;
typedef Vec3<double> Vec3d;

typedef std::vector<std::vector<ulong> > part_t;

double magic_polynomial(double t);



double bump(double t);
bool subset(double a,double b,double c,double d);
bool intersect(double a,double b,double c,double d);


Vec3d firstGaussDerivative(double s,Vec3d evalAt);
void secondGaussDerivative(double s, double theExp,Vec3d evalAt,
                           double* d2fdxx,
                           double* d2fdyy,
                           double* d2fdzz,
                           double* d2fdxy,
                           double* d2fdxz,
                           double* d2fdyz);



Matrix3d MakeMatrix3d(double a00, double a01, double a02,
                      double a10, double a11, double a12,
                      double a20, double a21, double a22);



#endif
