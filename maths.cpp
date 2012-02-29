

#include "maths.hpp"
#include "foreach.hpp"

#include <boost/array.hpp>

using namespace std;
using namespace boost;
using namespace SpinXML;


double sq(double x) {return x*x;}

array<array<ulong,3>,6> perm012 = {{
        {{0,1,2}},
        {{0,2,1}},
        {{1,0,2}},
        {{2,0,1}},
        {{1,2,0}},
        {{2,1,0}}
    }}; //Each shows up only once

array<array<ulong,3>,3> pem011 = {{
        {{0,1,1}},
        {{1,0,1}},
        {{1,1,0}}
    }}; //Each shows up 2! times

ulong pem000[1][3] = {
    {0,0,0}
}; //Each shows up 3! times

// 3! + 3*2! + 6 = (3^2)(2^2)/2
// 6 + 6 + 6 = 18  = 9*4/2 

template<typename T>
void applyPerm3(Vec3<T> in,const array<ulong,3> perm,Vec3<T>& out) {
    for(ulong i = 0; i < 3; i++) {
        out.set(perm[i],in[i]);
    }
}


part_t partitions_of_n(ulong n)  {
    part_t out;
    if(n == 0) {
        return part_t();
    }
    vector<ulong> trivial;
    trivial.push_back(n);
    out.push_back(trivial);


    //================================================================================//
    //Zoghbi and Stojmenovic's algorithm
    vector<ulong> x;
    x.resize(n);

    for(ulong i = 0; i < n; i++) {
        x[i] = 1;
    }
    x[0] = n;                
    long m  = 1;
    long h  = 1;

    while(x[0] != 1) {
        if(x[h-1] == 2) {
            m++;
            x[h-1] = 1;
            h--;
        } else {
            long r = x[h-1]-1;
            long t = m - h + 1;
            x[h-1] = r;
            while(t >= r) {
                h++;
                x[h-1] = r;
                t -= r;
            }
            if(t == 0) {
                m = h;
            } else {
                m = h + 1;
                if(t > 1) {
                    h++;
                    x[h-1] = t;
                }
            }
        }

        //Output
        vector<ulong> thisOut;
        thisOut.resize(m);
        for(long i = 0; i<m; i++) {
            thisOut[i] = x[i];
        }
        out.push_back(thisOut);

    }
    // End Zoghbi and Stojmenovic
    //================================================================================//

    return out;
}

ulong factorial(ulong n) {
    if(n == 0) return 1;
    if(n == 1) return 1;
    return n*factorial(n-1);
}

ulong triFactorial(ulong n, ulong m, ulong p) {
    return factorial(n)*factorial(m)*factorial(p);
}

//The number of ways of making 
ulong part_coef(vector<ulong> partition) {
    ulong coef = 1;
    ulong last = partition[0];
    ulong block_size = 0;

    for(ulong i = 0; i < partition.size(); i++) {
        block_size++;
        if(last != partition[i]) {
            last = partition[i];
            coef*=factorial(block_size);
            block_size = 1;
        }
        last = partition[i];
    }
    coef*=factorial(block_size);
    return coef;
}

GaussTaylorTerm::GaussTaylorTerm(ulong number)
    : termOrder(number) {
    part_t partitions = partitions_of_n(number);

    for(ulong n = 0; n <= termOrder; n++) {
        for(ulong m = 0; m < termOrder+1-n; m++) {
            ulong p = termOrder - m - n;

            Vec3ul nmp(n,m,p);

            MNL Pnmp = MNL::one();

            for(ulong j = 0; j < 3; j++) {
                for(ulong k = 0; k < nmp[j]; k++) {
                    Pnmp = hermite_recursion(Pnmp,j);
                }
            }

            mTerms.insert(pair<Vec3ul,MNL>(nmp,Pnmp));
            
        }
    }

}

MNL hermite_recursion(MNL mn,ulong dim) {
    return mn.take_derivative(dim) + mn.mult_by_power(dim,1)*(-2);
}


//Unit is the physical unit (dimensions of length) of nm. It
//should be equil to the total number of derivatives that have
//been taken
double eval_multinomial_with_dimension(MNL nm,long unit,double s,Vec3d xyz) {
    assert(s!= 0);
    assert(isfinite(xyz[0]));
    assert(isfinite(xyz[1]));
    assert(isfinite(xyz[2]));

    assert(isfinite(s));

    double total = 0;

    ulong size = nm.getSize();

    for(ulong i = 0; i < size; i++) {
        double term = 1;

        ulong ijk[3];
        long coef;

        nm.itor(i,&coef,ijk);
        if(coef == 0) {
            continue;
        }

        for(ulong j = 0; j<3; j++) {
            term*=pow(xyz[j],ijk[j]);
            assert(isfinite(term));
        }
        
        //If s and x,y,z must have the same units and so we can tell what
        //the power of s is in each term by chosing the power of s so that
        //the total power becomes unit.
        ulong term_power = ijk[0] + ijk[1] + ijk[2];

        long spower = unit - term_power;
        total += coef*term*pow(s,spower);
        assert(isfinite(total));        
    }
    assert(isfinite(total));
    return total;
}


double GaussTaylorTerm::eval(Vec3d r, Vec3d r0,double s) const {
    typedef pair<Vec3ul,MNL> t_thisPair;

    double dx = r[0] - r0[0];
    double dy = r[1] - r0[1];
    double dz = r[2] - r0[2];

    Vec3d rmr0(dx,dy,dz);

    double arg  = -sq(dx)- sq(dy)- sq(dz);
    double theExp = exp(arg/s*s);

    double total = 0;
    foreach(t_thisPair p,mTerms) {
        Vec3ul thisPart = p.first;
        MNL   u_nmp = p.second;

        long unit = -termOrder;
        //Is the derivative without the exponential part
        double poly   = eval_multinomial_with_dimension(u_nmp,unit,s,rmr0);

        assert(isfinite(poly));

        //cout << "facto " << pow(dx,thisPart[0])*pow(dy,thisPart[1])*pow(dz,thisPart[2]);
        //cout << "      " << poly;
        //cout << "      " << triFactorial(thisPart[0],thisPart[1],thisPart[2]) << endl;

        total += pow(dx,thisPart[0])*pow(dy,thisPart[1])*pow(dz,thisPart[2])
            *poly/triFactorial(thisPart[0],thisPart[1],thisPart[2]);
        assert(isfinite(total));
    }

    //Remember to multiply by the exponetial part
    double ret = total*theExp;
    assert(isfinite(ret));
    return ret;
}

/*
double GaussTaylorTerm::eval(Vec3d a, Vec3d x,double s) {
    typedef pair<Vec3ul,MNL> t_thisPair;

    double total = 0;

    foreach(t_thisPair p, mTerms) {
        Vec3l thisPart = p.first;
        MNL   u_nmp = p.second;

        if(thisPart[0] == thisPart[1] && thisPart[1] == thisPart[2]) {
            //They are all the same
            total+= u_nmp.eval(x) * 3* factorial(thisPart[0]);
        } else if (thisPart[0] != thisPart[1] && 
                   thisPart[1] != thisPart[2] && 
                   thisPart[0] != thisPart[2]) {
            //They are all different
            for(ulong i = 0; i < 6; i++) {
                Vec3d out(0,0,0);
                applyPerm3(x,perm012[i],out);
                total += u_nmp.eval(out);
            }
        } else {
            //Two are the same and one is different
            if(thisPart[0] == thisPart[1]) {
                //Shift the position of the part 2
                total += 2*u_nmp.eval(Vec3d(x[0],x[1],x[2]));
                total += 2*u_nmp.eval(Vec3d(x[0],x[2],x[1]));
                total += 2*u_nmp.eval(Vec3d(x[2],x[0],x[1]));
            } else if(thisPart[1] == thisPart[2]) {
                //Shift the position of the part 0                
                total += 2*u_nmp.eval(Vec3d(x[0],x[1],x[2]));
                total += 2*u_nmp.eval(Vec3d(x[1],x[0],x[2]));
                total += 2*u_nmp.eval(Vec3d(x[1],x[2],x[0]));
            } else {
                assert(thisPart[0] == thisPart[1]);
                //Shift the position of the part 1 
                total += 2*u_nmp.eval(Vec3d(x[1],x[0],x[2]));
                total += 2*u_nmp.eval(Vec3d(x[0],x[1],x[2]));
                total += 2*u_nmp.eval(Vec3d(x[0],x[2],x[1]));
            }
        }
    }
    return total;
    }*/


std::ostream& operator<< (std::ostream& out,const GaussTaylorTerm& term) {
    typedef pair<Vec3ul,MNL> t_thisPair;
    foreach(t_thisPair p, term.mTerms) {
        Vec3ul thisPart = p.first;
        MNL   u_nmp = p.second;

        out << noshowpos;
        out << "u(" << thisPart[0] << "," << thisPart[1] << "," << thisPart[2] << ") = " << u_nmp;
        out << showpos;

        out << endl;
    }
    return out;
}

vector<GaussTaylorTerm> terms;

const GaussTaylorTerm& getTerm(ulong i) {
    while(terms.size() < i+1) {
        terms.push_back(GaussTaylorTerm(terms.size()));
    }
    assert(i < terms.size());
    return terms[i];
}

double eval_gaussian_terms(Vec3d r, Vec3d r0,double s,ulong start) {
    double ret = 0;

    ulong i = start;
    double thisTerm;
    do {
        thisTerm = getTerm(i).eval(r,r0,s);
        assert(isfinite(thisTerm));

        ret += thisTerm;
        i++;
        cout << "ith term = " << thisTerm << " ith appox = " << ret << endl;
    } while(abs(thisTerm) > 1e-10);
    assert(isfinite(ret));
    return ret;
}

double gaussian_error_term_one(Vec3d r, Vec3d a, double s) {
    //Computes the exact remainder of the second order taylor
    //expansion of exp(-r2/s2) expanded around a

    double s2 = s*s;

    double x = r[0];               double x0 = a[0];
    double y = r[1];               double y0 = a[1];
    double z = r[2];               double z0 = a[2];
                                                   
    double x2 = x*x;               double x0_sq = x0*x0;
    double y2 = y*y;               double y0_sq = y0*y0;
    double z2 = z*z;               double z0_sq = z0*z0;

    double exp_x2 = exp(-x2/s2);   double exp_x0_sq = exp(-x0_sq/s2);  
    double exp_y2 = exp(-y2/s2);   double exp_y0_sq = exp(-y0_sq/s2);
    double exp_z2 = exp(-z2/s2);   double exp_z0_sq = exp(-z0_sq/s2);

    double int_x =  (exp_x2 + exp_x0_sq * (1 + 2 * x0*(x0-x)/s2) )   * exp(-(y0_sq + z0_sq)/s2);
    double int_y =  (exp_y2 + exp_y0_sq * (1 + 2 * y0*(y0-y)/s2) )   * exp(-(x0_sq + z0_sq)/s2);
    double int_z =  (exp_z2 + exp_z0_sq * (1 + 2 * z0*(z0-z)/s2) )   * exp(-(x0_sq + y0_sq)/s2);

    double int_x_y = (exp(-(x0_sq+y0_sq)/s2) - exp(-(x2+y0_sq)/s2) - exp(-(x0_sq+y2)/s2) + exp(-(x2+y2)/s2)) * exp_z0_sq;
    double int_y_z = (exp(-(y0_sq+z0_sq)/s2) - exp(-(y2+z0_sq)/s2) - exp(-(y0_sq+z2)/s2) + exp(-(y2+z2)/s2)) * exp_x0_sq;
    double int_x_z = (exp(-(x0_sq+z0_sq)/s2) - exp(-(x2+z0_sq)/s2) - exp(-(x0_sq+z2)/s2) + exp(-(x2+z2)/s2)) * exp_y0_sq;

    double int_x_y_z =
        -exp((x0_sq + y0_sq  + z0_sq) /s2)
        +exp((x0_sq + y0_sq  + z2)    /s2)
        -exp((x0_sq + y2     + z0_sq) /s2)
        +exp((x0_sq + y2     + z2)    /s2)
        -exp((x2    + y0_sq  + z0_sq) /s2)
        +exp((x2    + y0_sq  + z2)    /s2)
        -exp((x2    + y2     + z0_sq) /s2)
        +exp((x2    + y2     + z2)    /s2);
 
    return int_x + int_y + int_z + int_x_y + int_y_z + int_x_z + int_x_y_z;
}

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

Vec3d firstGaussDerivative(double s,Vec3d evalAt) {
        double s2 = s*s;
        double theExp = exp( -(evalAt.x()*evalAt.x() + 
                               evalAt.y()*evalAt.y() +
                               evalAt.z()*evalAt.z())/s2);

        //First derivative
        return Vec3d(-2*evalAt.x()/s2 * theExp,
                     -2*evalAt.y()/s2 * theExp,
                     -2*evalAt.z()/s2 * theExp);
}

void secondGaussDerivative(double s, double theExp,Vec3d evalAt,
                           double* d2fdxx,
                           double* d2fdyy,
                           double* d2fdzz,
                           double* d2fdxy,
                           double* d2fdxz,
                           double* d2fdyz) {
    double s2 = s*s;
    double s4 = s2*s2;


    *d2fdxx = (4*evalAt.x()*evalAt.x() - 2*s2)/s4 * theExp;
    *d2fdyy = (4*evalAt.y()*evalAt.y() - 2*s2)/s4 * theExp;
    *d2fdzz = (4*evalAt.z()*evalAt.z() - 2*s2)/s4 * theExp;

    *d2fdxy = (4*evalAt.x()*evalAt.y())/s4 * theExp;
    *d2fdxz = (4*evalAt.x()*evalAt.z())/s4 * theExp;
    *d2fdyz = (4*evalAt.y()*evalAt.z())/s4 * theExp;
}

/*
EulerAngles MatrixToEuler(const Matrix3d& rot) {
    Vector3d z_axis=Vector3d(0,0,1);
    Vector3d x_axis=Vector3d(1,0,0);
    z_axis=rot*z_axis;
    x_axis=rot*x_axis;
    
    double gamma=atan2(z_axis.y(),z_axis.x());
    double beta=atan2(sqrt(z_axis.x()*z_axis.x() + z_axis.y()*z_axis.y()),z_axis.z());

    //Use γ and β to rotate V2 back onto the point 0,0,1
    Quaterniond betaTwist(AngleAxisd(-beta, Vector3d(0,1,0)));
    Quaterniond gammaTwist (AngleAxisd(-gamma,Vector3d(0,0,1)));

    x_axis=(betaTwist*gammaTwist)*x_axis;

    double alpha = atan2(x_axis.y(),x_axis.x());
    if(alpha < 0 || alpha >= 2*PI)  alpha = alpha-2*PI*floor(alpha/(2*PI));
    if(alpha >= 2*PI) alpha = 0;
    if(beta < 0 || beta >= PI)    beta = beta-  PI*floor(beta/PI);
    if(gamma < 0 || gamma >= 2*PI)  gamma = gamma-2*PI*floor(gamma/(2*PI));
    if(gamma >= 2*PI) gamma = 0;
    return EulerAngles(alpha,beta,gamma);
}
*/

Matrix3d MakeMatrix3d(double a00, double a01, double a02,
                      double a10, double a11, double a12,
                      double a20, double a21, double a22) {
    Matrix3d m;
    m(0,0)=a00;             m(0,1)=a01;             m(0,2)=a02;
    m(1,0)=a10;             m(1,1)=a11;             m(1,2)=a12;
    m(2,0)=a20;             m(2,1)=a21;             m(2,2)=a22;

    return m;
}

AxRhomTensor tensorToAxRhom(Tensor tensor) {
    Matrix3d matrix = tensor.asMatrix();

    SpinXML::InteractionPayload pl = SpinXML::InteractionPayload(matrix);
    SpinXML::AxRhom      ar = pl.AsAxRhom();
    SpinXML::EulerAngles ea = ar.mOrient.GetAsEuler();

    AxRhomTensor axRhomTensor;
    axRhomTensor.ax    = ar.ax;
    axRhomTensor.rh    = ar.rh;
    axRhomTensor.alpha = ea.alpha;
    axRhomTensor.beta  = ea.beta;
    axRhomTensor.gamma = ea.gamma;

    return axRhomTensor;
}


Tensor axRhomToTensor(AxRhomTensor axRhomTensor) {
    double ax = axRhomTensor.ax;
    double rh = axRhomTensor.rh;

    double ev_xx = ( rh - ax/3)/2;
    double ev_yy = (-rh - ax/3)/2;
    double ev_zz = ax/3;

    SpinXML::EulerAngles ea(axRhomTensor.alpha,
                            axRhomTensor.beta,
                            axRhomTensor.gamma);

    Matrix3d dcm = ConvertToDCM(ea);
    
    Matrix3d eigenFrame = MakeMatrix3d(ev_xx,0     ,0,
                                       0,     ev_yy,0,
                                       0,     0     ,ev_zz);

    Matrix3d chi_tensor = dcm * eigenFrame * dcm.transpose();

    double chi_xx = chi_tensor(0,0);
    double chi_yy = chi_tensor(1,1);
    double chi_zz = chi_tensor(2,2);

    Tensor tensor;
    tensor.chi_1 = (2*chi_zz - chi_xx - chi_yy)/3.0;
    tensor.chi_2 = (chi_xx - chi_yy)/6;

    tensor.chi_xy = chi_tensor(0,1);
    tensor.chi_xz = chi_tensor(0,2);
    tensor.chi_yz = chi_tensor(1,2);

    return tensor;
}


Matrix3d Tensor::asMatrix() const {
    double chi_xx = -1/2.0*chi_1 + 3*chi_2;
    double chi_yy = -1/2.0*chi_1 - 3*chi_2;
    double chi_zz = chi_1;

    return MakeMatrix3d(chi_xx,chi_xy,chi_xz,
                        chi_xy,chi_yy,chi_yz,
                        chi_xz,chi_yz,chi_zz);
}
