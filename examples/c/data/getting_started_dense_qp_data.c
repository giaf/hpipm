/***************
* dim
***************/
/* nv */
int nv = 2;
/* ne */
int ne = 1;
/* nb */
int nb = 2;
/* ng */
int ng = 0;
/* nsb */
int nsb = 0;
/* nsg */
int nsg = 0;
/* ns */
int ns = 0;
/***************
* qp
***************/
/* H */
static double HH[] = {4.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, };
double *H = HH;
/* A */
static double AA[] = {1.000000000000000e+00, 1.000000000000000e+00, };
double *A = AA;
/* C */
static double CC[] = {};
double *C = CC;
/* idxb */
static int iidxb[] = {0, 1, };
int *idxb = iidxb;
/* g */
static double gg[] = {1.000000000000000e+00, 1.000000000000000e+00, };
double *g = gg;
/* zl */
static double zzl[] = {};
double *zl = zzl;
/* zu */
static double zzu[] = {};
double *zu = zzu;
/* b */
static double bb[] = {1.000000000000000e+00, };
double *b = bb;
/* lb */
static double llb[] = {-1.000000000000000e+00, -1.000000000000000e+00, };
double *lb = llb;
/* lb_mask */
static double llb_mask[] = {1.000000000000000e+00, 1.000000000000000e+00, };
double *lb_mask = llb_mask;
/* ub */
static double uub[] = {1.500000000000000e+00, 5.000000000000000e-01, };
double *ub = uub;
/* ub_mask */
static double uub_mask[] = {1.000000000000000e+00, 1.000000000000000e+00, };
double *ub_mask = uub_mask;
/* lg */
static double llg[] = {};
double *lg = llg;
/* lg_mask */
static double llg_mask[] = {};
double *lg_mask = llg_mask;
/* ug */
static double uug[] = {};
double *ug = uug;
/* ug_mask */
static double uug_mask[] = {};
double *ug_mask = uug_mask;
/* lls */
static double llls[] = {};
double *lls = llls;
/* lls_mask */
static double llls_mask[] = {};
double *lls_mask = llls_mask;
/* lus */
static double llus[] = {};
double *lus = llus;
/* lus_mask */
static double llus_mask[] = {};
double *lus_mask = llus_mask;
/* Zl */
static double ZZl[] = {};
double *Zl = ZZl;
/* Zu */
static double ZZu[] = {};
double *Zu = ZZu;
/* idxs_rev */
static int iidxs_rev[] = {-1, -1, };
int *idxs_rev = iidxs_rev;
/***************
* arg
***************/
/* mode */
int mode = 1;
/* iter_max */
int iter_max = 20;
/* alpha_min */
double alpha_min = 1.000000000000000e-12;
/* mu0 */
double mu0 = 1.000000000000000e+01;
/* tol_stat */
double tol_stat = 1.000000000000000e-12;
/* tol_eq */
double tol_eq = 1.000000000000000e-12;
/* tol_ineq */
double tol_ineq = 1.000000000000000e-12;
/* tol_comp */
double tol_comp = 1.000000000000000e-12;
/* reg_prim */
double reg_prim = 1.000000000000000e-15;
/* reg_dual */
double reg_dual = 1.000000000000000e-15;
/* warm_start */
int warm_start = 0;
/* pred_corr */
int pred_corr = 1;
/* split_step */
int split_step = 1;
