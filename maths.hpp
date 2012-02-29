
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
#include <spinxml/interaction.hpp>

//Of course, you shouldn't use namespaces in header files, but we'll
//probably be alright in this case.
using namespace Eigen;

typedef unsigned long ulong;

struct AxRhomTensor {
    AxRhomTensor() {}
    AxRhomTensor(double _ax,double _rh,double _alpha,double _beta, double _gamma)
        : ax(_ax),rh(_rh),alpha(_alpha),beta(_beta),gamma(_gamma) {
    }
    double ax;
    double rh;

    double alpha;
    double beta;
    double gamma;
};

struct Tensor {
    Tensor() {}
    Tensor(double _chi_1,double _chi_2,double _chi_xy,double _chi_xz, double _chi_yz)
        : chi_1(_chi_1),chi_2(_chi_2),chi_xy(_chi_xy),chi_xz(_chi_xz),chi_yz(_chi_yz) {
    }
    Eigen::Matrix3d asMatrix() const;
    double chi_1;
    double chi_2;

    double chi_xy;
    double chi_xz;
    double chi_yz;
};



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

template<class T> struct Multinomial;

template<class T>
std::ostream& operator <<(std::ostream& out,const Multinomial<T>& mn) {
    bool first = true;
    for(ulong i = 0; i < mn.getDim(2); i++) {
        for(ulong j = 0; j < mn.getDim(1); j++) {
            for(ulong k = 0; k < mn.getDim(0); k++) {
                ulong buff[3];
                buff[0] = k;
                buff[1] = j;
                buff[2] = i;
                T val = mn(buff);
                
                if(val == 0) {
                    continue;
                }
                if(first) {
                    out << " + ";
                    first = false;
                }
                out << val;
                if(k > 0) {
                    out << "x";
                    if(k > 1) {
                        out << "^" << k;
                    }
                }
                if(j > 0) {
                    out << "y";
                    if(j > 1) {
                        out << "^" << j;
                    }
                }
                if(i > 0) {
                    out << "z";
                    if(i > 1) {
                        out << "^" << i;
                    }
                }

            }
        }
    }

    return out;
}


template<class T>
struct Multinomial {
public:
    Multinomial(const ulong _size[3])
        : mData(new T[_size[0]*_size[1]*_size[2]]) {
        assert(_size[0]*_size[1]*_size[2] > 0);
        memcpy(sizes, _size,3*sizeof(ulong));
        for(ulong i = 0;i<_size[0]*_size[1]*_size[2]; i++) {
            assert(validI(i));
            mData[i] = 0;
        }
    }
    Multinomial(const Multinomial& other) 
        : mData(other.mData) {
        memcpy(sizes, other.sizes,3*sizeof(ulong));
    }
    Multinomial& operator=(const Multinomial& other) {
        mData = other.mData;
        memcpy(sizes, other.sizes,3*sizeof(ulong));
        return *this;
    }
    ~Multinomial() {
    }

    Multinomial copy() const {
        Multinomial cp(sizes);
        for(ulong i = 0; i < sizes[0]*sizes[1]*sizes[2]; i++) {
            assert(validI(i));
            assert(cp.validI(i));
            cp.mData[i] = mData[i];
        }
        return cp;
    }
    ulong getDim(ulong dim) const {
        return sizes[dim];
    }

    ulong getSize() const {
        return sizes[0]*sizes[1]*sizes[2];
    }

    void itor(ulong i, T* coef, ulong* ijk) const {
        assert(i < getSize());
        assert(coef != NULL);
        
        Multinomial::i2ijk(i,sizes,ijk);
        *coef = mData[i];
    }

    T& at(const ulong* ijk) const {
        return (*this)(ijk);
    }
    T& operator()(const ulong* ijk) const {
        for(ulong i = 0; i< 3; i++) {
            assert(ijk[i] >= 0);
            assert(ijk[i] <= sizes[i]);
        }

        return mData[Multinomial::ijk2i(sizes,ijk)];
    }
    double eval(Vec3d xyz) {
        assert(std::isfinite(xyz[0]));
        assert(std::isfinite(xyz[1]));
        assert(std::isfinite(xyz[2]));

        double total = 0;

        for(ulong i = 0; i<getSize(); i++) {
            ulong ijk[3];
            Multinomial::i2ijk(i,sizes,ijk);

            if(mData[i] == 0) {
                continue;
            }

            double term = mData[i];
            for(ulong j = 0; j<3; j++) {
                term*=pow(xyz[j],ijk[j]);
            }
            total += term;
        }
        assert(std::isfinite(total));
        return total;
    }

    Multinomial take_derivative(ulong dim) const {
        ulong newSizes[3];
        memcpy(newSizes,sizes,3*sizeof(ulong));
        newSizes[dim]--;
        if(newSizes[dim] == 0) {
            //The derivative of a constant is zero
            newSizes[0] = 1;
            newSizes[1] = 1;
            newSizes[2] = 1;
            return Multinomial(newSizes);
        }
        Multinomial ret(newSizes);
        ulong newSize = newSizes[0]*newSizes[1]*newSizes[2];

        for(ulong i=0; i < newSize; i++) {
            ulong ijk[3];
            i2ijk(i,newSizes,ijk);

            ulong sourceijk[3];
            memcpy(sourceijk,ijk,3*sizeof(ulong));
            
            sourceijk[dim]++;

            ulong source_i = ijk2i(sizes,sourceijk);

            assert(ret.validI(i));
            assert(validI(source_i));

            ret.mData[i] += mData[source_i] * sourceijk[dim];
        }
        return ret;
    }
    
    Multinomial operator+(const Multinomial& other) const {
        ulong maxSizes[3];
        for(ulong i = 0; i<3 ;i++) {
            maxSizes[i] = sizes[i] > other.sizes[i] ? sizes[i] : other.sizes[i];
        }
        Multinomial ret(maxSizes);

        for(ulong i=0; i < maxSizes[0]*maxSizes[1]*maxSizes[2]; i++) {
            ulong ijk[3];
            i2ijk(i,maxSizes,ijk);

            if(ijk[0] < sizes[0] && ijk[1] < sizes[1] && ijk[2] < sizes[2]) {
                ulong i2 = ijk2i(sizes,ijk);

                assert(ret.validI(i));
                assert(validI(i2));
                ret.mData[i] += mData[i2];
            }

            if(ijk[0] < other.sizes[0] && ijk[1] < other.sizes[1] && ijk[2] < other.sizes[2]) {
                ulong i2 = ijk2i(other.sizes,ijk);
                assert(ret.validI(i));
                assert(other.validI(i2));
                ret.mData[i] += other.mData[i2];
            }
        }
        return ret;
    }

    Multinomial operator*(T s) const {
        Multinomial ret(sizes);
        for(ulong i=0; i < sizes[0]*sizes[1]*sizes[2]; i++) {
            ret.mData[i] = s*mData[i];
        }
        return ret;
    }

    static Multinomial one() {
        ulong buff[3];

        buff[0] = buff[1] = buff[2] = 1;
        Multinomial ret(buff);

        ret.mData[0] = 1;
        return ret;
    }


    Multinomial mult_by_power(ulong dim,ulong power) const {
        ulong retSizes[3];
        memcpy(retSizes,sizes,sizeof(ulong)*3);

        retSizes[dim]+=power;
        Multinomial ret(retSizes);

        for(ulong i=0; i < sizes[0]*sizes[1]*sizes[2]; i++) {
            ulong ijk[3];
            i2ijk(i,sizes,ijk);

            ijk[dim] += power;

            ret.mData[ijk2i(retSizes,ijk)] = mData[i];
        }
        return ret;
    }

    friend std::ostream& operator<< <> (std::ostream& out,const Multinomial<T>& mn);
private:
    bool validI(ulong i) const {
        return i < sizes[0]*sizes[1]*sizes[2];
    }
    static ulong ijk2i(const ulong _sizes[3], const ulong ijk[3]) {
        return ijk[0] + ijk[1]*(_sizes[0]) + ijk[2]*_sizes[0]*_sizes[1];
    }
    static void i2ijk(ulong i,const ulong _sizes[3],ulong ijk[3]) {
        ijk[0] = i                         % _sizes[0];
        ijk[1] = (i/_sizes[0])             % _sizes[1];
        ijk[2] = (i/(_sizes[0]*_sizes[1])) % _sizes[2];
    }

    //Storage order: x dimension is the least significant, the z
    //dimension is the most significant
    ulong sizes[3];
    boost::shared_array<T> mData;
};

typedef Multinomial<long>   MNL;

class GaussTaylorTerm;

std::ostream& operator<< (std::ostream& out,const GaussTaylorTerm& term);

class GaussTaylorTerm {
public:
    GaussTaylorTerm(ulong number);

    //Evauate at x expanded around a
    double eval(Vec3d r, Vec3d r0,double s) const;

    friend std::ostream& operator<< (std::ostream& out,const GaussTaylorTerm& term);
private:
    GaussTaylorTerm();
    
    ulong termOrder;
    std::map<Vec3ul,MNL> mTerms;
};

typedef std::vector<std::vector<ulong> > part_t;

double eval_gaussian_terms(Vec3d a, Vec3d x,double s,ulong start);

double gaussian_error_term_one(Vec3d r, Vec3d a, double s);

MNL hermite_recursion(MNL mn,ulong dim);
part_t partitions_of_n(ulong n);
ulong factorial(ulong n);
ulong part_coef(std::vector<ulong> partition);

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


struct GaussDerivatives {
public:
    GaussDerivatives(double s,Vec3d xyz) {
        double sx = xyz.x();
        double sy = xyz.y();
        double sz = xyz.z();

        double s2 = s*s;
        double s4 = s2*s2;
        double s6 = s4*s2;

        double rho = exp(-(sx*sx + sy*sy + sz*sz)/s*s);
        //First derivative
        d_rho0_x = -2*sx/s2 * rho;
        d_rho0_y = -2*sy/s2 * rho;
        d_rho0_z = -2*sz/s2 * rho;

        //Second Derivative
        d2_rho0_xx = (4*sx*sx - 2*s2)/s4 * rho;
        d2_rho0_yy = (4*sy*sy - 2*s2)/s4 * rho;
        d2_rho0_zz = (4*sz*sz - 2*s2)/s4 * rho;

        d2_rho0_xy = (4*sx*sy)/s4 * rho;
        d2_rho0_xz = (4*sx*sz)/s4 * rho;
        d2_rho0_yz = (4*sy*sz)/s4 * rho;

        //Third Derivative
        d3_rho0_xxx = (12*sx*s2 - 8*sx*sx*sx)/s6 * rho;
        d3_rho0_yyy = (12*sy*s2 - 8*sy*sy*sy)/s6 * rho;
        d3_rho0_zzz = (12*sz*s2 - 8*sz*sz*sz)/s6 * rho;
               
        d3_rho0_xxy = (4*sy*(s2-2*sx*sx))/s6 * rho;
        d3_rho0_xxz = (4*sz*(s2-2*sx*sx))/s6 * rho;
               
        d3_rho0_yyx = (4*sx*(s2-2*sy*sy))/s6 * rho;
        d3_rho0_yyz = (4*sz*(s2-2*sy*sy))/s6 * rho;
               
        d3_rho0_zzx = (4*sx*(s2-2*sz*sz))/s6 * rho;
        d3_rho0_zzy = (4*sy*(s2-2*sz*sz))/s6 * rho;
               
        d3_rho0_xyz = -(8*sx*sy*sz)/s6 * rho;
    }
    //First derivative
    double d_rho0_x;
    double d_rho0_y;
    double d_rho0_z;

    //Second Derivative
    double d2_rho0_xx;
    double d2_rho0_yy;
    double d2_rho0_zz;

    double d2_rho0_xy;
    double d2_rho0_xz;
    double d2_rho0_yz;

    //Third Derivative
    double d3_rho0_xxx;
    double d3_rho0_yyy;
    double d3_rho0_zzz;
                      
    double d3_rho0_xxy;
    double d3_rho0_xxz;
                      
    double d3_rho0_yyx;
    double d3_rho0_yyz;
                      
    double d3_rho0_zzx;
    double d3_rho0_zzy;

    double d3_rho0_xyz;
};


Matrix3d MakeMatrix3d(double a00, double a01, double a02,
                      double a10, double a11, double a12,
                      double a20, double a21, double a22);


//-----------------------------------------------------------//
//------------------ Tensor Rotations -----------------------//
//-----------------------------------------------------------//

AxRhomTensor tensorToAxRhom(Tensor tensor);
Tensor axRhomToTensor(AxRhomTensor axRhomTensor);


#endif
