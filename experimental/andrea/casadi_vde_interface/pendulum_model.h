#ifndef EXAMPLES_C_PENDULUM_MODEL_PENDULUM_MODEL_H_
#define EXAMPLES_C_PENDULUM_MODEL_PENDULUM_MODEL_H_

int vdeFun(const double** arg, double** res, int* iw, double* w, int mem);

void VDE_fun_pendulum(const double* in, double* out,
                      int (*vde)(const double**, double**, int*, double*, int));

#endif  // EXAMPLES_C_PENDULUM_MODEL_PENDULUM_MODEL_H_
