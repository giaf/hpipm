#include "model.h"

int main(){

    printf("Testing CasADi-generated VDE evaluation:\n\n");
    const int NX = 4;
    const int NU = 1;

    // test input
    double in[NX + NU + NX*NX + NX*NU];
    double out[NX + NX*NX + NX*NU];

    for (int i = 0; i < NX + NU + NX*NX + NX*NU; i++) in[i] = i;

    printf("input:\n\n");
    for (int i = 0; i < NX;  i++) printf("x[%i] = %f\n", i, in[i]);
    printf("\n");
    for (int i = 0; i < NX*NX;  i++) printf("dx/dx[%i] = %f\n", i, in[i+NX]);
    printf("\n");
    for (int i = 0; i < NX*NU;  i++) printf("dx/du[%i] = %f\n", i, in[i+NX*NX]);
    printf("\n");
    for (int i = 0; i < NU;  i++) printf("u[%i] = %f\n", i, in[i+NX*NX+NX*NU]);
    printf("\n");

    VDE_fun_model(in, out, vdeFun);

    printf("output:\n\n");
    for (int i = 0; i < NX;  i++) printf("x[%i] = %f\n", i, out[i]);
    printf("\n");
    for (int i = 0; i < NX*NX;  i++) printf("dx/dx[%i] = %f\n", i, out[i+NX]);
    printf("\n");
    for (int i = 0; i < NX*NU;  i++) printf("dx/du[%i] = %f\n", i, out[i+NX*NX]);
    printf("\n");

    printf("-> done.\n");
}
