
#include <boost/function.hpp>
#include <vector>

/********************************************************************************
 * minimiser class
 ********************************************************************************/

template<typename T>
class Minimiser {
public:
	typedef std::vector<double> PList;
	Minimiser(boost::function<T(const PList&)> unpack,
			  boost::function<PList(const T&)> pack,  
			  boost::function<double(T)> min_funcion,		  
			  boost::function<void(const T&)> logger)		  
		: mUnpack(unpack),mPack(pack),mFMin(min_funcion),mLogger(logger) {
	}

	double minimise(T startingModel) {
		//Unpack the starting model
		PList plist = mPack(startingModel);
		unsigned long nParams = plist.size();

		//Initalise the minimiser
		gsl_multimin_fminimizer* gslmin =
			gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,nParams);

		gsl_multimin_function minfunc;
		minfunc.f = &Minimiser::minf_static;
		minfunc.n = nParams;
		minfunc.params = (void*)this;

		gsl_vector* vec = pList2GSLVec(plist);
		gsl_vector* step_size = gsl_vector_alloc(plist.size());
		
		for(unsigned long i = 0;i<plist.size();i++) {
			gsl_vector_set(step_size,i,plist[i] == 0 ? 1 : plist[i]*0.5);
		}

		gsl_multimin_fminimizer_set (gslmin, &minfunc, vec, step_size);
		for(unsigned long i = 0; i<1000;i++) {
			gsl_multimin_fminimizer_iterate (gslmin);
		}

        double min = gsl_multimin_fminimizer_minimum(gslmin);
		gsl_multimin_fminimizer_free(gslmin);
		gsl_vector_free(step_size);
        return min;
	}
private:
	double _minf(const gsl_vector * v) {
		T params = mUnpack(gSLVec2pList(v));
        double result = mFMin(params);
		mLogger(params);
		return result;
	}
	static double minf_static(const gsl_vector * v, void* _this) {
		return ((Minimiser*)_this)->_minf(v);
	}

	static gsl_vector* pList2GSLVec(PList pList) {
		gsl_vector* vec = gsl_vector_alloc (pList.size());
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

	boost::function<T(const PList&)> mUnpack;
	boost::function<PList(const T&)> mPack;
	boost::function<double(T)> mFMin;
	boost::function<void(T)> mLogger;	  
};
