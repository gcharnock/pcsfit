
#include "pointdev.hpp"
#include "model2.hpp"

#include <boost/shared_array.hpp>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cstring>

using namespace std;
using namespace boost;


typedef unsigned long ulong;

template<class T> struct Multinomial;

template<class T>
ostream& operator <<(ostream& out,const Multinomial<T>& mn) {
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
    double eval(const double xyz[3]) {
        double total = 0;

        for(ulong i = 0; i<sizes[0]*sizes[1]*sizes[2]; i++) {
            ulong ijk[3];
            Multinomial::i2ijk(i,sizes,ijk);

            double term = mData[i];
            for(ulong i = 0; i<3; i++) {
                term*=pow(xyz[i],ijk[i]);
            }
            total += term;
        }
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

    Multinomial u_plus_one_recusion(ulong dim,ulong n) const {
        Multinomial uprime = this->take_derivative(dim);

        Multinomial uprime_x2y2z2 =
            (uprime.mult_by_power(0,2) + 
             uprime.mult_by_power(1,2) + 
             uprime.mult_by_power(2,2));

        Multinomial twonaddfive = 
            this->mult_by_power(dim,1)*(-(2*n+5));

        return uprime_x2y2z2 + twonaddfive;
    }

    friend ostream& operator<< <> (ostream& out,const Multinomial<T>& mn);
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
    shared_array<T> mData;
};



//Evauates the numerator and derivatives of the point model
/*
  Derivatives come out in this order
  {0,0,0},
  
  {1,0,0},
  {0,1,0},
  {0,0,1},

  {1,1,0},
  {1,0,1},
  {0,1,1},

  {2,0,0},
  {0,2,0},
  {0,0,2}
*/

void eval_point_numerator(Vector3 evalAt,const double* pm,double out[9]) {
    double x = evalAt.x - pm[PARAM_X];
    double y = evalAt.y - pm[PARAM_Y];
    double z = evalAt.z - pm[PARAM_Z];

    double chi_1 =  pm[PARAM_CHI1]; 
    double chi_2 =  pm[PARAM_CHI2];
    double chi_xy = pm[PARAM_CHIXY];
    double chi_xz = pm[PARAM_CHIXZ];
    double chi_yz = pm[PARAM_CHIYZ];

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;

    double xy = x*y;
    double xz = x*z;
    double yz = y*z;

    double r2 = x2+y2+z2;



    out[0] = (r2-3*x2)*chi_1 + (z2-y2)*chi_2 + 6*(xy*chi_xy + xz*chi_xz + yz*chi_yz);

    out[1] = -4*x*chi_1 +             6*(y*chi_xy + z*chi_xz);
    out[2] =  2*y*chi_1 - 2*y*chi_2 + 6*(x*chi_xy + z*chi_yz);
    out[3] =  2*z*chi_1 + 2*z*chi_2 + 6*(x*chi_xz + y*chi_yz);

    out[4] = 6*chi_xy;
    out[5] = 6*chi_yz;
    out[6] = 6*chi_yz;

    out[7] = -4*chi_1;   
    out[8] =  2*chi_1 - 2*chi_2;
    out[9] =  2*chi_1 + 2*chi_2;

    for(unsigned long i = 0;i<9;i++) {
        assert(isfinite(out[i]));
    }
}

typedef Multinomial<long> MNL;


double u_nmp_over_v(ulong n,ulong m,ulong p, double x,double y,double z) {
    MNL u = MNL::one();
    for(ulong i = 0; i<n; i++) {
        u = u.u_plus_one_recusion(0,i);
    }
    for(ulong i = 0; i<m; i++) {
        u = u.u_plus_one_recusion(1,i+n);
    }
    for(ulong i = 0; i<p; i++) {
        u = u.u_plus_one_recusion(2,i+n+m);
    }

    ulong nTot = n + m + p;

    double vn = pow(x*x+y*y+z*z, (2.0*nTot+5.0)/2.0);

    double xyz[3];
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    return u.eval(xyz)/vn;
}

double eval_point_model_dev_xyz(const double* params,ulong ndx, ulong ndy, ulong ndz,Vector3 evalAt) {
    double x = evalAt.x - params[PARAM_X];
    double y = evalAt.y - params[PARAM_Y];
    double z = evalAt.z - params[PARAM_Z];

    double numerator_devs[9];
    eval_point_numerator(evalAt,params,numerator_devs);

    //There are ten possible partitions of 2, so we will just write
    //them out explicitly

    static const ulong two_partitions[10][3] = {
        {0,0,0},

        {1,0,0},
        {0,1,0},
        {0,0,1},

        {1,1,0},
        {1,0,1},
        {0,1,1},

        {2,0,0},
        {0,2,0},
        {0,0,2}
    };

    //Each of these partitions has an associated multinomial
    //coefficent NB: This must be in the same order as the
    //two_partitions array
    ulong multinomial_coef[10];

    multinomial_coef[0] = 1;

    multinomial_coef[1] = ndx;
    multinomial_coef[2] = ndy;
    multinomial_coef[3] = ndz;

    multinomial_coef[4] = ndx*ndy;
    multinomial_coef[5] = ndx*ndz;
    multinomial_coef[6] = ndy*ndz;

    multinomial_coef[7] = ndx*(ndx-1)/2;
    multinomial_coef[8] = ndy*(ndy-1)/2;
    multinomial_coef[9] = ndz*(ndz-1)/2;

    
    //Compute the derivative
    double total = 0;

    for(ulong i = 0;i<10;i++) {
        ulong ntwid = ndx - two_partitions[i][0];
        ulong mtwid = ndy - two_partitions[i][1];
        ulong ptwid = ndz - two_partitions[i][2];

        if(multinomial_coef[i] == 0) {
            continue;
        }
        total += multinomial_coef[i]*u_nmp_over_v(ntwid,mtwid,ptwid,x,y,z)*numerator_devs[i];
    }

    return total;
}

void multinomial_run_tests() {
    MNL u = MNL::one();
    MNL u200 = MNL::one();
    cout << "u000 = " << u << endl;

    for(ulong i = 1; i < 10; i++) {
        u = u.u_plus_one_recusion(0,i-1);
        cout << "u" << i << "00 = " << u << endl;

        if(i == 2) {
            u200 = u.copy();
        }
    }
    
    cout << "================================================================================" << endl;
    
    MNL u2 = u200;
    for(ulong i = 0; i < 5; i++) {
        u2 = u2.u_plus_one_recusion(1,i+2);
        cout << "u2" << (i+1) << "0 = " << u2 << endl;
    }

    cout << "================================================================================" << endl;
    double params[8];

    //Use integers for testing values because they work better with
    //mathematica
    long metal_x = 0,metal_y = 0, metal_z = 0,
        chi_1 = 1, chi_2 = 1, chi_xy = 1, chi_xz = 1, chi_yz = 1;

    long x = 4;
    long y = 5;
    long z = 4;

    params[PARAM_X]  = metal_x;
    params[PARAM_Y]  = metal_y;
    params[PARAM_Z]  = metal_z;

    params[PARAM_CHI1 ] = chi_1;
    params[PARAM_CHI2 ] = chi_2;
    params[PARAM_CHIXY] = chi_xy;
    params[PARAM_CHIXZ] = chi_xz;
    params[PARAM_CHIYZ] = chi_yz;

    double derivative = eval_point_model_dev_xyz(params,2,3,2,Vector3(4,5,4));
    cout << "derivative = " << derivative << endl;

    cout << "Expression to past into mathematica/wolfram alpha" << endl;
    cout << "D[((y^2+z^2-2*x^2)*" << chi_1
         << "+(z^2-y^2)*"         << chi_2
         << "+6*(x*y*" << chi_xy
         <<" + x*z*"   << chi_xz
         <<" + y*z*"   << chi_yz
         << "))*(x^2+y^2+z^2)^(-5/2),{x,2},{y,3},{z,2}] /."
         << " { x->" << (x-metal_x)
         << ",y-> " << (y-metal_y)
         << ", z-> " << (z-metal_z) << "}" << endl;
}

