

#include "tests.hpp"
#include "model2.hpp"
#include "pointdev.hpp"
#include "fit.hpp"
#include "maths.hpp"
#include "foreach.hpp"

#include <boost/bind.hpp>
#include <vector>
#include <iostream>

using namespace SpinXML;
using namespace std;

RandomDist dist;

void check_minimum(PRNG& prng,const Model* model,Multithreader<fdf_t>* pool,const ModelOptions* modelOptions) {
    cout << "Evautating the error function for a perfect match (should be zero)" << endl;

    for(unsigned long i = 0;i<10;i++) {
        double* params = (double*)alloca(model->size*sizeof(double));
        for(unsigned long j=0;j<model->size;j++) {
            params[j] = dist(prng);
        }
        Dataset dataset;
        random_data(prng,*model,params,5,&dataset,modelOptions);
        
        ErrorContext context;
        context.dataset = &dataset;
        context.params  = params;
        context.model   = model;
        context.pool    = pool;
        context.rescale = false;

        double  error,n_error;
        double* gradient = (double*)alloca(model->size*sizeof(double));
        double* n_gradient = (double*)alloca(model->size*sizeof(double));

        eval_error(&context,params,&error,gradient);
        numerical_error_derivative(&context,&n_error,n_gradient);

        cout << "run " << i << ", error = " << error << ", n_error = " << n_error << " (should be zero)" << endl;
        cout << "Analytic: "; for(unsigned long j = 0;j < model->size;j++) {cout << gradient[j] << " ";} cout << endl;
        cout << "Numeric: "; for(unsigned long j = 0;j < model->size;j++) {cout << n_gradient[j] << " ";}  cout << endl;
    }
}

void check_derivative (PRNG& prng,const Model* model,const ModelOptions* modelOptions) {
    cout << "Checking the analytic and numerical derivaties match on " << model->name << endl;
    
    //Chekcs the analytic and numerical derivates are equil
    for(unsigned long i = 0;i<3;i++) {
        double result,point_result;
        double* params             = (double*)alloca(model->size*sizeof(double));
        double* gradient           = (double*)alloca(model->size*sizeof(double));
		double* numerical_gradient = (double*)alloca(model->size*sizeof(double));

        for(unsigned long j = 0; j< model->size; j++) {
            params[j] = dist(prng);
        }

		Vec3d evalAt(dist(prng),dist(prng),dist(prng));

        model->modelf(evalAt,params,&result,gradient,modelOptions);
        point_model.modelf(evalAt,params,&point_result,NULL,modelOptions);

        numerical_derivative(evalAt,model,params,numerical_gradient,modelOptions);

        double x = evalAt.x() - params[PARAM_X];
        double y = evalAt.y() - params[PARAM_Y];
        double z = evalAt.z() - params[PARAM_Z];

        cout << "evalAt = (" << evalAt.x() << "," << evalAt.y() << "," << evalAt.z() << ") r = "
             << sqrt(x*x+y*y+z*z) << endl;
        cout << "Result = " << result << " (point result = " << point_result <<  ")" << endl;
        for(unsigned long i=0;i<model->size;i++) {
			cout << name_param(i) << " = " << params[i] <<
				" grad = " << gradient[i]
				 << " ceneral Differenceing gradient of " << numerical_gradient[i] << endl;
        }
        cout << "================================================================================" << endl;
    }
}

void check_error_derivate(PRNG prng,const Model* model,Multithreader<fdf_t>* pool,const ModelOptions* modelOptions) {
    cout << "Checking the analytic and numerical derivaties match for the error functional" << endl;

    double* params             = (double*)alloca(model->size*sizeof(double));
    double* params2            = (double*)alloca(model->size*sizeof(double));
    double* gradient           = (double*)alloca(model->size*sizeof(double));
    double* numerical_gradient = (double*)alloca(model->size*sizeof(double));

    //Chekcs the analytic and numerical derivates are equil
    for(unsigned long i = 1;i<4;i++) {
        double result;

		//We need to generate two models, or else the gradient of the
		//error function will always be zero (and the numerical
		//gradient will be tiny and based off numerical errors)
        for(unsigned long j = 0; j< model->size; j++) {
            params[j] = dist(prng);
            params2[j] = dist(prng);
        }

        Dataset dataset;

		cout << "Creating a random molecule of " << 5*i << " spins" << endl;
        random_data(prng,point_model,params,5*i,&dataset,modelOptions);

        ErrorContext context;
        context.dataset = &dataset;
        context.model   = model;
        context.params  = params2;
        context.pool    = pool;

        eval_error(&context,params2,&result,gradient);

        numerical_error_derivative(&context,&result,numerical_gradient);

        cout << "Result = " << result << endl;
        for(unsigned long i=0;i<model->size;i++) {
			cout << name_param(i) << " = " << params[i];

			cout.setf(ios::floatfield,ios::scientific);
			cout <<	" grad = " << gradient[i]
				 << " ceneral Differenceing gradient of " << numerical_gradient[i] << endl;
        }
        cout << "================================================================================" << endl;
    }
}

bool do_convergence_on_iterate(const ErrorContext* context,unsigned long i,gsl_multimin_fdfminimizer* minimizer) {
	double fx = gsl_multimin_fdfminimizer_minimum(minimizer);
	gsl_vector* g = gsl_multimin_fdfminimizer_gradient(minimizer);
    gsl_vector* x = gsl_multimin_fdfminimizer_x(minimizer);

	double norm = 0;
	for(unsigned long j = 0; j < g->size; j++) {
		norm += gsl_vector_get(g,j)*gsl_vector_get(g,j);
	}

    cout << i << "\t";
    for(unsigned long j = 0;j< x->size;j++) {cout << gsl_vector_get(x,j) << "\t";}
	cout << "\tf(x) = " << fx << "\t|grad| = " << norm << endl;

	return true;
}

void do_convergence(PRNG& prng,const Model* model,Multithreader<fdf_t>* pool,const ModelOptions* modelOptions) {
	unsigned long size = model->size;

    for(unsigned long i = 1;i<4;i++) {
        double* params_real        = (double*)alloca(size*sizeof(double));
        double* params_start       = (double*)alloca(size*sizeof(double));
        double* params_opt         = (double*)alloca(size*sizeof(double));

        for(unsigned long j=0;j<size;j++) {
            params_real[j]  = dist(prng);
            params_start[j] = params_real[j] + 0.1*dist(prng);
        }

        Dataset dataset;

        cout << "Generating a distom molecule with " << 10*i << " spins" << endl;
        random_data(prng,point_model,params_real,5*i,&dataset,modelOptions);

        cout << "Real params:" << endl;
        for(unsigned long j = 0;j<model->size;j++) {cout << params_real[j] << " ";}
        cout << endl;

        ErrorContext context;
        context.dataset = &dataset;
        context.params  = params_start;
        context.model   = model;
        context.pool    = pool;

		double errorStart = 666;
        double errorFinal = 666;
		double* gradient  = (double*)alloca(size*sizeof(double));

		eval_error(&context,params_start,&errorStart,gradient);

        do_fit_with_grad(&context,params_opt,&errorFinal,do_convergence_on_iterate);
        for(unsigned long j = 0;j < 8; j++) {
            cout << name_param(j) << ": real = " << params_real[j] << " start = " << params_start[j]
                 << " final = " << params_opt[j] << endl;
        }
		cout.setf(ios::floatfield,ios::scientific);

        cout << "The inital error was " << errorStart << " and the final error was: " << errorFinal << endl;

        cout << "================================================================================" << endl;
    }
}

void test_gaussian(PRNG& prng,const ModelOptions* modelOptions) {
	RandomDist rand;

    for(unsigned long i = 1; i<4; i++) {
        //Generate a random dataset using the Gaussian and Gaussian
        //test models. We expect them to match

        double* params = (double*)alloca(gaussian_model.size*sizeof(double));

        for(unsigned long j=0;j<gaussian_model.size;j++) {
            params[j]  = dist(prng);
        }
        unsigned long natoms = i*10;

        Nuclei nuclei;
        Vals valsG;
        Vals valsTest;

        nuclei.resize(natoms);
        valsG.resize(natoms);
        valsTest.resize(natoms);

        for(unsigned long j = 0; j < natoms;j++) {
            Vec3d pos = Vec3d(rand(prng),rand(prng),rand(prng));
            nuclei[j] = pos;
            eval_gaussian(        pos,params,&(valsG[j])   ,NULL,modelOptions);
            eval_gaussian_testing(pos,params,&(valsTest[j]),NULL,modelOptions);
        }


        for(unsigned long j = 0; j<natoms; j++) {
            double x = nuclei[j].x();
            double y = nuclei[j].y();
            double z = nuclei[j].z();

            double stddev = params[PARAM_STDDEV];
            double r_singularity = sqrt(x*x+y*y+z*z);

            cout << "(" << (r_singularity/stddev) << "," << valsG[j] << "," << valsTest[j] << ")" << endl;
        }
        cout << endl;
        cout << "================================================================================" << endl;
    }
}

//Perform a sanity check on the models. Given teh same paramiter set
//the gaussian model should converged to the point model as stddev
//tends to 0
void testModel(PRNG& prng,Multithreader<fdf_t>* pool,const ModelOptions* modelOptions) {
    cout << "================================================================================" << endl;

    //Check the gaussian model;
    //test_gaussian(prng);

    //check_derivative (prng,&point_model,modelOptions);
    //PRNG prng_copy = prng;
    check_derivative (prng,&gaussian_model,modelOptions);
    //check_derivative (prng_copy,&gaussian_model_num_dev,modelOptions);

    //cout << "Evaulating the analytic and numerical derivatives of the error functional" << endl;
    //check_error_derivate(prng,&point_model   ,pool);
    //check_error_derivate(prng,&gaussian_model,pool);

    //check_minimum(prng,&point_model   ,pool);
    //check_minimum(prng,&gaussian_model,pool);
    
	//do_convergence(prng,&point_model   ,pool);
	//do_convergence(prng,&gaussian_model,pool);
}

//================================================================================//
void multinomial_run_tests(const ModelOptions* modelOptions) {
    MNL u = MNL::one();
    MNL u200 = MNL::one();
    cout << "u000 = " << u << endl;

    for(ulong i = 1; i < 10; i++) {
        u = u_plus_one_recusion(u,0,i-1);
        cout << "u" << i << "00 = " << u << endl;

        if(i == 2) {
            u200 = u.copy();
        }
    }
    
    cout << "================================================================================" << endl;
    
    MNL u2 = u200;
    for(ulong i = 0; i < 5; i++) {
        u2 = u_plus_one_recusion(u2,1,i+2);
        cout << "u2" << (i+1) << "0 = " << u2 << endl;
    }

    cout << "================================================================================" << endl;
    
    MNL u3 = MNL::one();
    for(ulong i = 0; i < 5; i++) {
        u3 = u_plus_one_recusion(u3,2,i);
        cout << "u00" << (i+1) << " = " << u3 << endl;
    }


    cout << "================================================================================" << endl;
    double params[POINT_SIZE];

    //Use integers for testing values because they work better with
    //mathematica
    long metal_x = 0,metal_y = 0, metal_z = 0,
        chi_1 = 1, chi_2 = 2, chi_xy = -2, chi_xz = 8, chi_yz = -9;

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

    double derivative = eval_point_model_dev_xyz(params,2,3,2,Vec3d(4,5,4));
    cout << "derivative = " << derivative << endl;

    cout << "Expression to past into mathematica/wolfram alpha" << endl;
    cout << "D[((y^2+z^2-2*x^2)*" << chi_1
         << "+(z^2-y^2)*"         << chi_2
         << "+6*(x*y*" << chi_xy
         <<" + x*z*"   << chi_xz
         <<" + y*z*"   << chi_yz
         << "))*(x^2+y^2+z^2)^(-5/2),{x,2},{y,3},{z,2}] /."
         << " { x->" << (x-metal_x)
         << ",y->"   << (y-metal_y)
         << ", z-> " << (z-metal_z) << "}" << endl;
    cout << "================================================================================" << endl;
    Vec3d evalAt = Vec3d(4,4,4);

    //Evauate the point model using the well tested method, to compare the terms
    double value;
    double* gradient = (double*)alloca(POINT_SIZE*sizeof(double));
    eval_point(evalAt,params,&value,gradient,modelOptions);

    cout << "First term (evauated at (4,5,4) )" << endl;
    cout << eval_point_model_dev_xyz(params,0,0,0,evalAt) << endl;
    cout << "Evauated by eval_point = " << value << endl;
    cout << endl;

    cout << "Checking derivatives match eval_point. They should have oposite signs" << endl;

    cout << eval_point_model_dev_xyz(params,1,0,0,evalAt) << " =? " <<  gradient[PARAM_X] << endl;
    cout << eval_point_model_dev_xyz(params,0,1,0,evalAt) << " =? " <<  gradient[PARAM_Y] << endl;
    cout << eval_point_model_dev_xyz(params,0,0,1,evalAt) << " =? " <<  gradient[PARAM_Z] << endl;
    cout << endl;

    cout << "Second derivatives (Total should be zero)" << endl;
    double sigma_xx = eval_point_model_dev_xyz(params,2,0,0,evalAt);
    cout << "sigma_xx = " << sigma_xx << endl;
    double sigma_yy = eval_point_model_dev_xyz(params,0,2,0,evalAt);
    cout << "sigma_yy = " << sigma_yy << endl;
    double sigma_zz = eval_point_model_dev_xyz(params,0,0,2,evalAt);
    cout << "sigma_zz = " << sigma_zz << endl;
    cout << "Total = " << sigma_xx + sigma_yy + sigma_zz << endl;
}


double testPoly6(Vec3d r) {
    //6th order term of the taylor expansion of exp(-x2-y2-z2)
    double x = r[0];
    double y = r[1];
    double z = r[2];

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;

    double x4 = x2*x2;
    double y4 = y2*y2;
    double z4 = z2*z2;

    double x6 = x4*x2;
    double y6 = y4*y2;
    double z6 = z4*z2;
    return -1/6.0*x6 - 1/6.0*y6 - 1/2.0*y4*z2 - 1/2.0*y2*z4 - 1/6.0*z6
        - 1/2.0*(y2 +z2)*x4 - 1/2.0*(y4 + 2*y2*z2 + z4)*x2;
}


void testMaths(PRNG& prng) {
	RandomDist rand;

    //multinomial_run_tests();

    cout << "================================================================================" << endl;
    cout << "Checking we can convert tensors from ax-rhom-euler to our custom format and back again." << endl;
    
    for(ulong i = 0; i < 10; i++) {
        cout << "================================================================================" << endl;

        AxRhomTensor axRhomTensor;

        double ax = axRhomTensor.ax =  rand(prng);  //Axiality is always positive
        double rh = axRhomTensor.rh = -(axRhomTensor.ax)*rand(prng); //Rhombicity is always negative

        double alpha = axRhomTensor.alpha = rand(prng) * 2 * M_PI;
        double beta  = axRhomTensor.beta  = rand(prng) * M_PI;
        double gamma = axRhomTensor.gamma = rand(prng) * 2 * M_PI;
        
        cout << "Inital Random Tensor" << endl;
        cout << "ax = " << ax << endl;
        cout << "rh = " << rh << endl;

        cout << "alpha = " << alpha/(2*M_PI)*360 << endl;
        cout << "beta = "  << beta /(2*M_PI)*360 << endl;
        cout << "gamma = " << gamma/(2*M_PI)*360 << endl;
        cout << endl;



        {
            cout << "DCM matrix constructed directly from Euler Angles" << endl;
            Matrix3d rotAlpha = MakeMatrix3d(cos(alpha),-sin(alpha),0,
                                             sin(alpha), cos(alpha),0,
                                             0         , 0         ,1);

            //ZYZ
            /*Matrix3d rotBeta  = MakeMatrix3d( cos(beta),0, sin(beta),
                                              0,        1, 0,
                                              -sin(beta),0, cos(beta));*/
            //ZXZ
            Matrix3d rotBeta  = MakeMatrix3d(1, 0,        0,
                                             0, cos(beta),-sin(beta),
                                             0,sin(beta),cos(beta));



            Matrix3d rotGamma = MakeMatrix3d(cos(gamma),-sin(gamma),0,
                                             sin(gamma), cos(gamma),0,
                                             0         , 0         ,1);

            Matrix3d dcm = rotGamma*rotBeta*rotAlpha;

            cout << dcm << endl;
            cout << endl;

            cout << "DCM Matrix converted back to euler angles" << endl;
            EulerAngles ea = ConvertToEuler(dcm);
            cout << ea.alpha*180/M_PI << ", " << ea.beta*180/M_PI << ", " << ea.gamma*180/M_PI << endl;
            cout << endl;
        }

        Tensor tensor = axRhomToTensor(axRhomTensor);

        energy e_xx = - ax/6.0 + rh/2.0;
        energy e_yy = - ax/6.0 - rh/2.0;
        energy e_zz =   ax/3.0;
        cout << "Eigenvalue - Euler Angle Representation" << endl;
        cout << "xx = " << e_xx << endl;
        cout << "yy = " << e_yy << endl;
        cout << "zz = " << e_zz << endl;
        cout << endl;


        cout << "Intermediate Representation Tesnor:" << endl;

        cout << "chi_1 = " << tensor.chi_1 << endl;
        cout << "chi_2 = " << tensor.chi_2 << endl;

        cout << "chi_xy = " << tensor.chi_xy << endl;
        cout << "chi_xz = " << tensor.chi_xz << endl;
        cout << "chi_yz = " << tensor.chi_yz << endl;
        cout << endl;

        Matrix3d matrix = tensor.asMatrix();

        AxRhomTensor axRhomTensor2 = tensorToAxRhom(tensor);

        cout << "Matrix Tensor Representation" << endl;
        cout << matrix << endl;
        cout << endl;


        cout << "Eigenvalue - dcm representation" << endl;

        EigenSolver<Matrix3d> solver(matrix, true); //true => compute eigenvectors
        Vector3d v1,v2,v3;//Eigenvectors

        double e1=real(solver.eigenvalues()[0]);
        double e2=real(solver.eigenvalues()[1]);
        double e3=real(solver.eigenvalues()[2]);

        Matrix3d eigenVectors = solver.eigenvectors().real();
        Matrix3d so3 = eigenVectors * eigenVectors.determinant();

        cout << e1 << ", " << e2 << ", " << e3 << endl;
        cout << so3 << endl;
        cout << "determinant = " << so3.determinant() << endl;
        cout << endl;

        cout << "Sorted Eigenvalue - dcm representation" << endl;

        double xx = e1,yy = e2,zz = e3;
        double tmp;
        
        Orientation o(so3);
        if(zz < yy) {
            tmp=zz; zz=yy; yy=tmp;
            //To get back to the original orientation, rotate by one right
            //angle in the positive direction about x and invert about x-y
            //plane. But becaues the tensor is an axis aligned ellipse,
            //mirroring about a coordinate plane is the identity.

            Orientation rotation(AngleAxisd(M_PI/2,Vector3d(1,0,0)));
            o = o * rotation;
        }
     	if(zz < xx) {
            tmp=zz; zz=xx; xx=tmp;
            //Rotate about y in the positive direction

            Orientation rotation(AngleAxisd(M_PI/2,Vector3d(0,1,0)));
            o = o * rotation;
        }
     	if(yy < xx) {
            tmp=xx; xx=yy; yy=tmp;
            //Rotate about z in the postive direction.

            Orientation rotation(AngleAxisd(M_PI/2,Vector3d(0,0,1)));
            o = o * rotation;
        }
        Matrix3d sortedSo3 = o.GetAsMatrix();
        cout << xx << ", " << yy << ", " << zz << endl;
        cout << sortedSo3 << endl;
        cout << "determinant = " << sortedSo3.determinant() << endl;
        cout << endl;

        cout << "SortedSo3 matrix used to reconstruct the interaction matrix" << endl;
        cout << sortedSo3*MakeMatrix3d(e1,0 ,0 ,
                                       0 ,e2,0 ,
                                       0 ,0 ,e3) * sortedSo3.transpose() << endl;
        cout << endl;

        cout << "Sorted so3 matrix as euler angles" << endl;
        EulerAngles ea = ConvertToEuler(sortedSo3);
        cout << "alpha = " << ea.alpha*180/M_PI << endl;
        cout << "beta  = " << ea.beta *180/M_PI << endl;
        cout << "gamma = " << ea.gamma*180/M_PI << endl;
        cout << endl;

        cout << "Symmetary Group" << endl;
        vector<Matrix3d> s2s2s2;
        s2s2s2.push_back(MakeMatrix3d(1,0,0,
                                      0,1,0,
                                      0,0,1));

        s2s2s2.push_back(MakeMatrix3d(-1,0,0,
                                      0,1,0,
                                      0,0,1));

        s2s2s2.push_back(MakeMatrix3d(1,0,0,
                                      0,-1,0,
                                      0,0,1));

        s2s2s2.push_back(MakeMatrix3d(1,0,0,
                                      0,1,0,
                                      0,0,-1));

        s2s2s2.push_back(MakeMatrix3d(-1,0,0,
                                      0,-1,0,
                                      0,0,1));

        s2s2s2.push_back(MakeMatrix3d(1,0,0,
                                      0,-1,0,
                                      0,0,-1));

        s2s2s2.push_back(MakeMatrix3d(-1,0,0,
                                      0,1,0,
                                      0,0,-1));

        s2s2s2.push_back(MakeMatrix3d(-1,0,0,
                                      0,-1,0,
                                      0,0,-1));

        for(unsigned long i = 0; i< s2s2s2.size(); i++) {
            Matrix3d m = sortedSo3*s2s2s2[i];
            if(m.determinant() < 0 ) {
                continue;
            }
            
            EulerAngles ea = ConvertToEuler(m);
            cout << "  alpha = " << ea.alpha*180/M_PI;
            cout << " beta  = " << ea.beta *180/M_PI;
            cout << " gamma = " << ea.gamma*180/M_PI << endl;
            cout << endl;

            
        }



        cout << "Final Random Tensor" << endl;
        cout << "ax = " << axRhomTensor2.ax << endl;
        cout << "rh = " << axRhomTensor2.rh << endl;

        cout << "alpha = " << axRhomTensor2.alpha/(2*M_PI)*360 << endl;
        cout << "beta = "  << axRhomTensor2.beta /(2*M_PI)*360 << endl;
        cout << "gamma = " << axRhomTensor2.gamma/(2*M_PI)*360 << endl;
        cout << endl;

        cout << "--------" << endl;
    }

    return;

    cout << "================================================================================" << endl;
    for(ulong i = 1; i < 8; i++) {
        part_t part_n = partitions_of_n(i);

        cout << "Paritions of " << i << " are: " << endl;
        foreach(const vector<ulong>& thisPart,part_n) {
            foreach(ulong n,thisPart) {
                cout << n << " ";
            }
            cout << endl;
        }
        cout << "----------------------------------------" << endl;
    }
    cout << "================================================================================" << endl;

    for(ulong i = 1; i < 8; i++) {
        cout << i << "! = " << factorial(i) << " ";
    }
    cout << endl;
    

    cout << "================================================================================" << endl;

    for(ulong i = 0; i < 5; i++) {
        GaussTaylorTerm term(i);
        cout << term << endl;
        cout << "================================================================================" << endl;        
    }

    cout << "================================================================================" << endl;
    {
        Vec3d r0 = Vec3d(1,1,1);
        Vec3d r = Vec3d(0,0,0);


        for(ulong i = 0; i < 7; i++) {
            GaussTaylorTerm term(i);
            cout << "Term " << i << " = " << term.eval(r,r0,2.0) << endl;
            cout << "================================================================================" << endl;
        }
    }
    cout << "Taylor approximating exp at random locations about random points" << endl;
    for(ulong i = 0; i < 10; i++) {
        double x = rand(prng);
        double y = rand(prng);
        double z = rand(prng);

        double x0 = x + 0.01*rand(prng);
        double y0 = y + 0.01*rand(prng);
        double z0 = z + 0.01*rand(prng);

        double dx = x - x0;
        double dy = y - y0;
        double dz = z - z0;

        double s = rand(prng);

        Vec3d r = Vec3d(x,y,z);
        Vec3d r0 = Vec3d(x0,y0,z0);

        cout << eval_gaussian_terms(r,r0,s,0) << " | " << exp(-(dx*dx + dy*dy + dz*dz)/(s*s)) << endl;;
    }
    
    cout << "================================================================================" << endl;

    cout << "Testing exact evauation of exp(r) - exp(r0) - (\\partial exp)(r0) via lagrange remainders" << endl;
    cout << "Lagrange Remainder | Via difference" << endl;

    cout << "Random Input" << endl;

    for(ulong i = 0; i < 10; i++) {
        //There are three important points to consider here. The
        //center of the gaussian, the point we are talyor expanding
        //around and the point we want to evaluate at. We will always
        //center the distribution on the origin, expand around r0 and
        //evaluate at r.

        //Generate input data
        double x = rand(prng);
        double y = rand(prng);
        double z = rand(prng);

        double s = rand(prng);

        double x0 = x + 0.1*rand(prng);
        double y0 = y + 0.1*rand(prng);
        double z0 = z + 0.1*rand(prng);

        //Calcuate some useful sub-expressions
        
        double r_sq = x*x + y*y + z*z;
        double r0_sq = x0*x0 + y0*y0 + z0*z0;

        //Calcuate the exponetial and its derivatives

        double a_coef = 1/(s*s);

        double norm = pow(M_PI*s*s,-1.5);
        double theExp  = exp(-a_coef*r_sq);
        double theExp0 = exp(-a_coef*r0_sq);

        double f  = norm*theExp;
        double f0 = norm*theExp0;


        double common_factor = -norm*2*theExp0/s;
        double drhoBydx_dx   = common_factor * x0 * (x-x0);
        double drhoBydy_dy   = common_factor * y0 * (y-y0);
        double drhoBydz_dz   = common_factor * z0 * (z-z0);

        Vec3d r  = Vec3d(x, y, z);
        Vec3d r0 = Vec3d(x0,y0,z0);


        cout << f << endl << f0 << endl;
        cout << drhoBydx_dx << " " << drhoBydy_dy << " " << drhoBydz_dz << endl;

        cout << "r/s = (" << x << "," << y << "," << z
             << ") | r0/s = (" << x0 << "," << y0 << "," << z0 << ");" << endl;
        cout << "r/s=" << sqrt(r_sq)/s <<  " r0_sq/s = " << sqrt(r0_sq)/s << " | " 
             << norm*gaussian_error_term_one(r,r0,s) << " | " 
             << (f - f0 - drhoBydx_dx - drhoBydy_dy - drhoBydz_dz)  << endl;
        cout << endl;
    }

    cout << endl;
    cout << "Known input" << endl;

    cout << "gaussian_error_term_one(Vec3d(1,2,3),Vec3d(-4,6,2),2) = " 
         << gaussian_error_term_one(Vec3d(1,2,3),Vec3d(-4,6,2),2) << endl;
    
}
