
#include <vector>
#include "data.hpp"
#include "maths.hpp"

double eval_point_model_dev_xyz(const double* params,ulong ndx, ulong ndy, ulong xdz,Vec3d evalAt);


template<typename T>
Multinomial<T> u_plus_one_recusion(Multinomial<T> mn, ulong dim, ulong n) {
    Multinomial<T> uprime = mn.take_derivative(dim);

    Multinomial<T> uprime_x2y2z2 =
        uprime.mult_by_power(0,2) + 
        uprime.mult_by_power(1,2) + 
        uprime.mult_by_power(2,2);

    Multinomial<T> twonaddfive = 
        mn.mult_by_power(dim,1)*(-(2*n+5));

    return uprime_x2y2z2 + twonaddfive;
}
