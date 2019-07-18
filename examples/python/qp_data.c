/***************
* dim
***************/
/* N */
int N = 5;
/* nx */
static int nnx[] = {2, 2, 2, 2, 2, 2, };
int *nx = nnx;
/* nu */
static int nnu[] = {1, 1, 1, 1, 1, 0, };
int *nu = nnu;
/* nbx */
static int nnbx[] = {2, 0, 0, 0, 0, 2, };
int *nbx = nnbx;
/* nbu */
static int nnbu[] = {0, 0, 0, 0, 0, 0, };
int *nbu = nnbu;
/* ng */
static int nng[] = {0, 0, 0, 0, 0, 0, };
int *ng = nng;
/* nsbx */
static int nnsbx[] = {0, 0, 0, 0, 0, 0, };
int *nsbx = nnsbx;
/* nsbu */
static int nnsbu[] = {0, 0, 0, 0, 0, 0, };
int *nsbu = nnsbu;
/* nsg */
static int nnsg[] = {0, 0, 0, 0, 0, 0, };
int *nsg = nnsg;
/***************
* qp
***************/
/* A */
static double A0[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };
static double A1[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };
static double A2[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };
static double A3[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };
static double A4[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };
static double *AA[] = {A0, A1, A2, A3, A4, };
double **hA = AA;
/* B */
static double B0[] = {0.000000000000000e+00, 1.000000000000000e+00, };
static double B1[] = {0.000000000000000e+00, 1.000000000000000e+00, };
static double B2[] = {0.000000000000000e+00, 1.000000000000000e+00, };
static double B3[] = {0.000000000000000e+00, 1.000000000000000e+00, };
static double B4[] = {0.000000000000000e+00, 1.000000000000000e+00, };
static double *BB[] = {B0, B1, B2, B3, B4, };
double **hB = BB;
/* b */
static double b0[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double b1[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double b2[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double b3[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double b4[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double *bb[] = {b0, b1, b2, b3, b4, };
double **hb = bb;
/* Q */
static double Q0[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q1[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q2[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q3[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q4[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q5[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double *QQ[] = {Q0, Q1, Q2, Q3, Q4, Q5, };
double **hQ = QQ;
/* S */
static double S0[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S1[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S2[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S3[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S4[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S5[] = {};
static double *SS[] = {S0, S1, S2, S3, S4, S5, };
double **hS = SS;
/* R */
static double R0[] = {1.000000000000000e+00, };
static double R1[] = {1.000000000000000e+00, };
static double R2[] = {1.000000000000000e+00, };
static double R3[] = {1.000000000000000e+00, };
static double R4[] = {1.000000000000000e+00, };
static double R5[] = {};
static double *RR[] = {R0, R1, R2, R3, R4, R5, };
double **hR = RR;
/* r */
static double r0[] = {0.000000000000000e+00, };
static double r1[] = {0.000000000000000e+00, };
static double r2[] = {0.000000000000000e+00, };
static double r3[] = {0.000000000000000e+00, };
static double r4[] = {0.000000000000000e+00, };
static double r5[] = {};
static double *rr[] = {r0, r1, r2, r3, r4, r5, };
double **hr = rr;
/* q */
static double q0[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q1[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q2[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q3[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q4[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q5[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double *qq[] = {q0, q1, q2, q3, q4, q5, };
double **hq = qq;
/* idxbu */
static int idxbu0[] = {};
static int idxbu1[] = {};
static int idxbu2[] = {};
static int idxbu3[] = {};
static int idxbu4[] = {};
static int idxbu5[] = {};
static int *iidxbu[] = {idxbu0, idxbu1, idxbu2, idxbu3, idxbu4, idxbu5, };
int **hidxbu = iidxbu;
/* lbu */
static double lbu0[] = {};
static double lbu1[] = {};
static double lbu2[] = {};
static double lbu3[] = {};
static double lbu4[] = {};
static double lbu5[] = {};
static double *llbu[] = {lbu0, lbu1, lbu2, lbu3, lbu4, lbu5, };
double **hlbu = llbu;
/* ubu */
static double ubu0[] = {};
static double ubu1[] = {};
static double ubu2[] = {};
static double ubu3[] = {};
static double ubu4[] = {};
static double ubu5[] = {};
static double *uubu[] = {ubu0, ubu1, ubu2, ubu3, ubu4, ubu5, };
double **hubu = uubu;
/* idxbx */
static int idxbx0[] = {0, 1, };
static int idxbx1[] = {};
static int idxbx2[] = {};
static int idxbx3[] = {};
static int idxbx4[] = {};
static int idxbx5[] = {0, 1, };
static int *iidxbx[] = {idxbx0, idxbx1, idxbx2, idxbx3, idxbx4, idxbx5, };
int **hidxbx = iidxbx;
/* lbx */
static double lbx0[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double lbx1[] = {};
static double lbx2[] = {};
static double lbx3[] = {};
static double lbx4[] = {};
static double lbx5[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double *llbx[] = {lbx0, lbx1, lbx2, lbx3, lbx4, lbx5, };
double **hlbx = llbx;
/* ubx */
static double ubx0[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double ubx1[] = {};
static double ubx2[] = {};
static double ubx3[] = {};
static double ubx4[] = {};
static double ubx5[] = {-0.000000000000000e+00, -0.000000000000000e+00, };
static double *uubx[] = {ubx0, ubx1, ubx2, ubx3, ubx4, ubx5, };
double **hubx = uubx;
/* C */
static double C0[] = {};
static double C1[] = {};
static double C2[] = {};
static double C3[] = {};
static double C4[] = {};
static double C5[] = {};
static double *CC[] = {C0, C1, C2, C3, C4, C5, };
double **hC = CC;
/* D */
static double D0[] = {};
static double D1[] = {};
static double D2[] = {};
static double D3[] = {};
static double D4[] = {};
static double D5[] = {};
static double *DD[] = {D0, D1, D2, D3, D4, D5, };
double **hD = DD;
/* lg */
static double lg0[] = {};
static double lg1[] = {};
static double lg2[] = {};
static double lg3[] = {};
static double lg4[] = {};
static double lg5[] = {};
static double *llg[] = {lg0, lg1, lg2, lg3, lg4, lg5, };
double **hlg = llg;
/* ug */
static double ug0[] = {};
static double ug1[] = {};
static double ug2[] = {};
static double ug3[] = {};
static double ug4[] = {};
static double ug5[] = {};
static double *uug[] = {ug0, ug1, ug2, ug3, ug4, ug5, };
double **hug = uug;
/* Zl */
double **hZl;
/* Zu */
double **hZu;
/* zl */
double **hzl;
/* zu */
double **hzu;
/* idxs */
double **hidxs;
/* lls */
double **hlls;
/* lus */
double **hlus;
/***************
* arg
***************/
/* mode */
int mode = 1;
/* iter_max */
int iter_max = 30;
/* alpha_min */
double alpha_min = 1.000000000000000e-12;
/* mu0 */
double mu0 = 1.000000000000000e+04;
/* tol_stat */
double tol_stat = 1.000000000000000e-04;
/* tol_eq */
double tol_eq = 1.000000000000000e-05;
/* tol_ineq */
double tol_ineq = 1.000000000000000e-05;
/* tol_comp */
double tol_comp = 1.000000000000000e-05;
/* reg_prim */
double reg_prim = 1.000000000000000e-12;
/* warm_start */
int warm_start = 0;
/* pred_corr */
int pred_corr = 1;
/* ric_alg */
int ric_alg = 1;
