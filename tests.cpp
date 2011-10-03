
#include "tests.hpp"
#include "model.hpp"
#include "model2.hpp"

#include <iostream>

using namespace std;


//Perform a sanity check on the models. Given teh same paramiter set
//the gaussian model should converged to the point model as stddev
//tends to 0
void testModel(long seed) {
	PRNG prng(seed);  //We should actually be recycling prng rather
					  //than the seed
	RandomDist rand;

	PointModel pm = PointModel::randomModel(seed+10);
    double pm2[8];

    pm2[PARAM_X]     = pm.metal.x;
    pm2[PARAM_Y]     = pm.metal.y;
    pm2[PARAM_Z]     = pm.metal.z;
           
    pm2[PARAM_CHI1]  = rand(prng);
    pm2[PARAM_CHI2]  = rand(prng);
    pm2[PARAM_CHIXY] = rand(prng);
    pm2[PARAM_CHIXZ] = rand(prng);
    pm2[PARAM_CHIYZ] = rand(prng);



	cout << "================================================================================" << endl;

	cout << "Checking the gaussian is normalised" << endl;
	cout << "TODO" << endl;

	cout << "================================================================================" << endl;
	cout << "Point Model = " << pm << endl;

	GaussModel gm1;

	gm1.ax = pm.ax;
	gm1.rh = pm.rh;

	gm1.metal = pm.metal;
	gm1.setEulerAngles(pm.angle_x,pm.angle_y,pm.angle_z);

	gm1.stddev = 1;

	GaussModel gm05 = gm1;
	gm05.stddev = 0.5;

	GaussModel gm01 = gm1;
	gm01.stddev = 0.1;

	GaussModel gm005 = gm1;
	gm005.stddev = 0.05;

	GaussModel gm001 = gm1;
	gm001.stddev = 0.01;

	GaussModel gm0005 = gm1;
	gm0005.stddev = 0.005;

	GaussModel gm0001 = gm1;
	gm0001.stddev = 0.001;

	GaussModel gm00005 = gm1;
	gm00005.stddev = 0.0005;


	cout << "================================================================================" << endl;
	for(unsigned long i = 0;i<0;i++) {
		double x = rand(prng);
		double y = rand(prng);
		double z = rand(prng);

		double dx = pm.metal.x - x;
		double dy = pm.metal.y - y;
		double dz = pm.metal.z - z;

		double r = sqrt(dx*dx + dy*dy + dz*dz);

		cout << "Point selected = ("<< x << "," << y << "," << z << ")" 
			 << "Metal-point distance is " << r << endl;

		cout << "Point Model = " << pm.eval(x,y,z) << endl;
		cout << "gm1 Model = " << gm1.eval(x,y,z) << endl;
		cout << "gm05 Model = " << gm05.eval(x,y,z) << endl;
		cout << "gm01 Model = " << gm01.eval(x,y,z) << endl;
		cout << "gm005 Model = " << gm005.eval(x,y,z) << endl;
		cout << "gm001 Model = " << gm001.eval(x,y,z) << endl;
		cout << "gm0005 Model = " << gm0005.eval(x,y,z) << endl;
		cout << "gm0001 Model = " << gm0001.eval(x,y,z) << endl;
		cout << "gm00005 Model = " << gm00005.eval(x,y,z) << endl;
		cout << "================================================================================" << endl;
	}

    cout << "Point Model 2 = "; for(unsigned long i=0;i<8;i++){cout << pm2[i] << " ";} cout << endl;

    cout << "================================================================================" << endl;

    for(unsigned long i = 0;i<10;i++) {
		pm2[PARAM_X]     = rand(prng);
		pm2[PARAM_Y]     = rand(prng);
		pm2[PARAM_Z]     = rand(prng);
           
		pm2[PARAM_CHI1]  = rand(prng);
		pm2[PARAM_CHI2]  = rand(prng);
		pm2[PARAM_CHIXY] = rand(prng);
		pm2[PARAM_CHIXZ] = rand(prng);
		pm2[PARAM_CHIYZ] = rand(prng);

        double result;
        double gradient[8];
		double numerical_gradient[8];

		Vector3 evalAt(rand(prng),rand(prng),rand(prng));

        eval_point(evalAt,pm2,&result,gradient);

        numerical_derivative(evalAt,pm2,eval_point,8,numerical_gradient);

        cout << "Result = " << endl;
        for(unsigned long i=0;i<8;i++) {
			cout << name_param(POINT_PARAM(i)) << " = " << pm2[i] <<
				" grad = " << gradient[i]
				 << " ceneral Differenceing gradient of " << numerical_gradient[i] << endl;
        }
        cout << "================================================================================" << endl;
    }

    for(unsigned long i = 0;i<1;i++) {
        double gm2[9];
		gm2[PARAM_X]      = rand(prng);
		gm2[PARAM_Y]      = rand(prng);
		gm2[PARAM_Z]      = rand(prng);
           
		gm2[PARAM_CHI1]   = rand(prng);
		gm2[PARAM_CHI2]   = rand(prng);
		gm2[PARAM_CHIXY]  = rand(prng);
		gm2[PARAM_CHIXZ]  = rand(prng);
		gm2[PARAM_CHIYZ]  = rand(prng);
        
        gm2[PARAM_STDDEV] = abs(rand(prng));

		Vector3 evalAt(rand(prng),rand(prng),rand(prng));
        
        double result;
        double gradient[9];
        double numerical_gradient[9];

        eval_gaussian(evalAt,gm2,&result,gradient);

        numerical_derivative(evalAt,gm2,eval_gaussian,9,numerical_gradient);

        cout << "Result = " << result << endl;
        for(unsigned long i=0;i<8;i++) {
			cout << name_param(POINT_PARAM(i)) << " = " << gm2[i] <<
				" grad = " << gradient[i]
				 << " ceneral Differenceing gradient of " << numerical_gradient[i] << endl;
        }
        cout << "stddev = " << gm2[8] <<
            " grad = " << gradient[8]
             << " ceneral Differenceing gradient of " << numerical_gradient[8] << endl;

        cout << "================================================================================" << endl;
    }
}


