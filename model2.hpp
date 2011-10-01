
#include <string>

typedef void (*ModelF) (const double*,double*,double*);


enum POINT_PARAM {
    PARAM_X,    
    PARAM_Y,    
    PARAM_Z,    
                
    PARAM_CHI1, 
    PARAM_CHI2, 
    PARAM_CHIXY,
    PARAM_CHIXZ,
    PARAM_CHIYZ,
};

void eval_point(const double pm[8],double* value, double gradient[8]);
std::string name_param(POINT_PARAM param);

struct Model {
    
};
