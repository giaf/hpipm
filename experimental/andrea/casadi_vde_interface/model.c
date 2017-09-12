#include "model.h"

void VDE_fun_model(const double* in, double* out,
                      int (*vde)(const double**, double**, int*, double*,
                                 int)) {

    const int NX = 4;
    const int NU = 1;
    const double* x = in;
    const double* Su = in + NX;
    const double* Sx = in + NX + NU * NX;
    const double* u = in + NX + NX * (NX + NU);

    double* x_out = out;
    double* Su_out = out + NX;
    double* Sx_out = out + NX + NU * NX;

    const double* casadi_arg[4];
    double* casadi_res[3];

    casadi_arg[0] = x;
    casadi_arg[1] = Sx;
    casadi_arg[2] = Su;
    casadi_arg[3] = u;

    casadi_res[0] = x_out;
    casadi_res[1] = Sx_out;
    casadi_res[2] = Su_out;

    int* iw = 0;
    double* w = 0;
    int mem = 0;

    vde(casadi_arg, casadi_res, iw, w, mem);
}
