
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>

#include "../../include/hpipm_d_ocp_qp_dim.h"
#include "../../include/hpipm_d_ocp_qp.h"
#include "../../include/hpipm_d_ocp_qp_sol.h"
#include "../../include/hpipm_d_ocp_qp_ipm.h"

#define QP_HORIZON 5

int main() {

    double A[] = {1, 0, 1, 1};
    double B[] = {0, 1};
    double b[] = {0, 0};

    double Q[] = {1, 0, 0, 1};
    double S[] = {0, 0};
    double R[] = {1};
    double q[] = {1, 1};
    double r[] = {0};

    double x0[] = {1, 1};
    int idxb0[] = {1, 2};

    int nx[] = {2, 2, 2, 2, 2, 2};
    int nu[] = {1, 1, 1, 1, 1, 0};
    int nb[] = {2, 0, 0, 0, 0, 0};
    int ng[] = {0, 0, 0, 0, 0, 0};
    int ns[] = {0, 0, 0, 0, 0, 0};
    int nbx[] = {2, 0, 0, 0, 0, 0};
    int nbu[] = {0, 0, 0, 0, 0, 0};

    int N = QP_HORIZON;

    /************************************************
    * ocp qp dim
    ************************************************/

	int dim_size = d_memsize_ocp_qp_dim(N);
	printf("\ndim size = %d\n", dim_size);
	void *dim_mem = malloc(dim_size);

	struct d_ocp_qp_dim dim;
	d_create_ocp_qp_dim(N, &dim, dim_mem);
	d_cvt_int_to_ocp_qp_dim(N, nx, nu, nbx, nbu, ng, ns, &dim);

    double *hA[] = {A, A, A, A, A};
    double *hB[] = {B, B, B, B, B};
    double *hb[] = {b, b, b, b, b};
    double *hQ[] = {Q, Q, Q, Q, Q, Q};
    double *hS[] = {S, S, S, S, S, S};
    double *hR[] = {R, R, R, R, R};
    double *hq[] = {q, q, q, q, q, q};
    double *hr[] = {r, r, r, r, r, r};
    int *hidxb[] = {idxb0};
    double *hlb[] = {x0};
    double *hub[] = {x0};

    /************************************************
    * ocp qp
    ************************************************/

	int qp_size = d_memsize_ocp_qp(&dim);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

	struct d_ocp_qp qp;
	d_create_ocp_qp(&dim, &qp, qp_mem);

    d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, 
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &qp);

    /************************************************
    * ocp qp sol
    ************************************************/

	int qp_sol_size = d_memsize_ocp_qp_sol(&dim);
	printf("\nqp sol size = %d\n", qp_sol_size);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_ocp_qp_sol qp_sol;
	d_create_ocp_qp_sol(&dim, &qp_sol, qp_sol_mem);

    /************************************************
    * ipm arg
    ************************************************/

	int ipm_arg_size = d_memsize_ocp_qp_ipm_arg(&dim);
	printf("\nipm arg size = %d\n", ipm_arg_size);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_ocp_qp_ipm_arg arg;
	d_create_ocp_qp_ipm_arg(&dim, &arg, ipm_arg_mem);
	d_set_default_ocp_qp_ipm_arg(&arg);

	int ipm_size = d_memsize_ocp_qp_ipm(&dim, &arg);
	printf("\nipm size = %d\n", ipm_size);
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qp_ipm_workspace workspace;
	d_create_ocp_qp_ipm(&dim, &arg, &workspace, ipm_mem);

	int hpipm_return; // 0 normal; 1 max iter

	int rep, nrep=1000;

	struct timeval tv0, tv1;

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
    {
		hpipm_return = d_solve_ocp_qp_ipm(&qp, &qp_sol, &arg, &workspace);
    }

	gettimeofday(&tv1, NULL); // stop

}
