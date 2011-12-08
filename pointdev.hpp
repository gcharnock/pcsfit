
#include <vector>
#include "data.hpp"

void multinomial_run_tests();

double eval_point_model_dev_xyz(const double* params,ulong ndx, ulong ndy, ulong xdz,Vector3 evalAt);
