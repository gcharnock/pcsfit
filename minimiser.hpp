
#include <functional>
#include <vector>

/********************************************************************************
 * minimiser class
 * A thin OO wrapper around the gsl_multimin* functions and nothing more
 ********************************************************************************/

template<typename T>
class MinimiserBase {
public:
	typedef std::vector<double> PList;

    MinimiserBase(boost::function<double(const T&)> min_funcion,
                  const PList& small_step_sizes,
                  T startingModel)
        : mFMin(min_funcion) {

		//Unpack the starting model
		pList = T::pack(startingModel);
		nParams = pList.size();

        startingModelVec = gsl_vector_alloc(nParams);
        pList2GSLVec(pList,startingModelVec);

        step_size = gsl_vector_alloc(nParams);
        pList2GSLVec(small_step_sizes,step_size);        
    }
    virtual ~MinimiserBase() {
        gsl_vector_free(step_size);
        gsl_vector_free(startingModelVec);
    }
    
    virtual std::pair<double,T> iterate() = 0;
protected:
	static gsl_vector* pList2GSLVec(PList pList,gsl_vector* vec) {
		for(unsigned long i = 0; i<pList.size();i++) {
			gsl_vector_set(vec,i,pList[i]);
		}
		return vec;
	}
	static PList gSLVec2pList(const gsl_vector* vec) {
		PList retVal;
		for(size_t i = 0;i<vec->size;i++) {
			retVal.push_back(gsl_vector_get(vec,i));
		}
		return retVal;
	}

	boost::function<double(T)> mFMin;
    PList pList;
    unsigned long nParams;
    gsl_vector* startingModelVec;
    gsl_vector* step_size;
};

template<typename T>
class SimplexMinimiser : public MinimiserBase<T> {
public:
    typedef std::vector<double> PList;
	SimplexMinimiser(boost::function<double(const T&)> min_funcion,
                     const PList& small_step_sizes,
                     T startingModel)		  
		: MinimiserBase<T>(min_funcion,small_step_sizes,startingModel) {
        gslmin = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,this->nParams);

		minfunc.f = &SimplexMinimiser::minf_static;
		minfunc.n = this->nParams;
		minfunc.params = (void*)this;
		gsl_multimin_fminimizer_set (gslmin, &minfunc, this->startingModelVec, this->step_size);
	}

    virtual ~SimplexMinimiser() {
        gsl_multimin_fminimizer_free(gslmin);
    }

protected:
    virtual std::pair<double,T> iterate() {
        gsl_multimin_fminimizer_iterate (gslmin);

        double min = gsl_multimin_fminimizer_minimum(gslmin);
        gsl_vector* minVec = gsl_multimin_fminimizer_x(gslmin);
        T best_model = T::unpack(this->gSLVec2pList(minVec));

        return std::pair<double,T>(min,best_model);
	}
private:
	double f(const gsl_vector * v) {
		T params = T::unpack(this->gSLVec2pList(v));
        double result = mFMin(params);
        assert(std::isfinite(result));
		return result;
	}
	static double minf_static(const gsl_vector * v, void* _this) {
		return ((SimplexMinimiser*)_this)->f(v);
	}

    gsl_multimin_function minfunc;
    gsl_multimin_fminimizer* gslmin;
};


template<typename T>
class GradientMinimiser : public MinimiserBase<T> {
public:
	typedef std::vector<double> PList;
	typedef boost::function<double(const T&)> ObjectiveFunction;
	typedef boost::function<double(const T&,bool,PList&)> GradientFunction;
	GradientMinimiser(ObjectiveFunction min_funcion, GradientFunction withGrad,
					  const PList& small_step_sizes, T startingModel)		  
		: MinimiserBase<T>(min_funcion,small_step_sizes,startingModel),mWithGrad(withGrad) {
        //Setup the function to be minimised
		minfunc.f = &GradientMinimiser::f_static;
		minfunc.df = &GradientMinimiser::df_static;
		minfunc.fdf = &GradientMinimiser::fdf_static;
		minfunc.n = this->nParams;
		minfunc.params = (void*)this;

		//Initalise the minimiser
		//gsl_multimin_fdfminimizer_vector_bfgs2
		gslmin = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr,this->nParams);
        gsl_multimin_fdfminimizer_set(gslmin,&minfunc,this->startingModelVec,1,0.1);
	}

    virtual ~GradientMinimiser() {
		gsl_multimin_fdfminimizer_free(gslmin);
    }

    virtual std::pair<double,T> iterate() {
        gsl_multimin_fdfminimizer_iterate(gslmin);

        double min = gsl_multimin_fdfminimizer_minimum(gslmin);
        gsl_vector* minVec = gsl_multimin_fdfminimizer_x(gslmin);
        T best_model = T::unpack(this->gSLVec2pList(minVec));

        return std::pair<double,T>(min,best_model);
	}
private:
	double f(const gsl_vector * v) {
		T params = T::unpack(this->gSLVec2pList(v));
        double result = mFMin(params);
        assert(std::isfinite(result));
		return result;
	}
	void df(const gsl_vector *v, gsl_vector *df) {
        /*gsl_vector* vprime = gsl_vector_alloc(v->size);
        cout << "grad =";

        for(unsigned long i = 2;i < v->size;i++) {
            double df_by_di = 0;
            double h = this->step_size->data[i];

            gsl_vector_memcpy(vprime,v);

            vprime->data[i] = v->data[i] - 2*h;
            df_by_di += f(vprime);

            vprime->data[i] = v->data[i] - 1*h;
            df_by_di -= 8*f(vprime);

            vprime->data[i] = v->data[i] + 1*h;
            df_by_di += 8*f(vprime);

            vprime->data[i] = v->data[i] + 2*h;
            df_by_di -= f(vprime);

            cout << " " << df_by_di;

            df_by_di/=(12*h);
            df->data[i] = df_by_di;
            /*
			vprime->data[i] = v->data[i]+h;
            df_by_di += f(vprime);

            vprime->data[i] = v->data[i]-h;
            df_by_di -= f(vprime);

            df_by_di/=(2*h);
            cout << " " << df_by_di;
			*/
		//            df->data[i] = df_by_di;
		//}
        //cout << endl;
        //gsl_vector_free(vprime);

		std::vector<double> g; g.resize(this->nParams);
		T params = T::unpack(this->gSLVec2pList(v));
		mWithGrad(params,true,g);
		this->pList2GSLVec(g,df);
	}


	static double f_static(const gsl_vector * v, void* _this) {
		return ((GradientMinimiser*)_this)->f(v);
	}
	static void df_static(const gsl_vector *v, void* _this, gsl_vector *df) {
		((GradientMinimiser*)_this)->df(v,df);
	}
	static void fdf_static(const gsl_vector * v, void* _this,double *f, gsl_vector *df) {
        *f = ((GradientMinimiser*)_this)->f(v); 
        ((GradientMinimiser*)_this)->df(v, df);
	}

	GradientFunction mWithGrad;
    gsl_multimin_fdfminimizer* gslmin;
    gsl_multimin_function_fdf minfunc;
};
