int N = 5;
static int nnx[] = {2, 2, 2, 2, 2, 2};
static int nnu[] = {1, 1, 1, 1, 1, 0};
static int nnbx[] = {2, 2, 2, 2, 2, 2};
static int nnbu[] = {1, 1, 1, 1, 1, 0};
static int nng[] = {0, 0, 0, 0, 0, 2};
static int nns[] = {0, 0, 0, 0, 0, 2};

static double A0[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };
static double A1[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };
static double A2[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };
static double A3[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };
static double A4[] = {1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.000000000000000e+00, };

static double B0[] = {0.000000000000000e+00, 1.000000000000000e+00, };
static double B1[] = {0.000000000000000e+00, 1.000000000000000e+00, };
static double B2[] = {0.000000000000000e+00, 1.000000000000000e+00, };
static double B3[] = {0.000000000000000e+00, 1.000000000000000e+00, };
static double B4[] = {0.000000000000000e+00, 1.000000000000000e+00, };

static double b0[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double b1[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double b2[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double b3[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double b4[] = {0.000000000000000e+00, 0.000000000000000e+00, };

static double Q0[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q1[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q2[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q3[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q4[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };
static double Q5[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };

static double R0[] = {1.000000000000000e+00, };
static double R1[] = {1.000000000000000e+00, };
static double R2[] = {1.000000000000000e+00, };
static double R3[] = {1.000000000000000e+00, };
static double R4[] = {1.000000000000000e+00, };
static double R5[] = {};

static double S0[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S1[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S2[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S3[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S4[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double S5[] = {};

static double q0[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q1[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q2[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q3[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q4[] = {1.000000000000000e+00, 1.000000000000000e+00, };
static double q5[] = {1.000000000000000e+00, 1.000000000000000e+00, };

static double r0[] = {0.000000000000000e+00, };
static double r1[] = {0.000000000000000e+00, };
static double r2[] = {0.000000000000000e+00, };
static double r3[] = {0.000000000000000e+00, };
static double r4[] = {0.000000000000000e+00, };
static double r5[] = {};

static int idxb0[] = {0, 1, 2, };
static int idxb1[] = {0, 1, 2, };
static int idxb2[] = {0, 1, 2, };
static int idxb3[] = {0, 1, 2, };
static int idxb4[] = {0, 1, 2, };
static int idxb5[] = {0, 1, };

static double lb0[] = {-5.000000000000000e-01, 1.000000000000000e+00, 1.000000000000000e+00, };
static double lb1[] = {-5.000000000000000e-01, -5.000000000000000e+00, -5.000000000000000e+00, };
static double lb2[] = {-5.000000000000000e-01, -5.000000000000000e+00, -5.000000000000000e+00, };
static double lb3[] = {-5.000000000000000e-01, -5.000000000000000e+00, -5.000000000000000e+00, };
static double lb4[] = {-5.000000000000000e-01, -5.000000000000000e+00, -5.000000000000000e+00, };
static double lb5[] = {-5.000000000000000e+00, -5.000000000000000e+00, };

static double ub0[] = {5.000000000000000e-01, 1.000000000000000e+00, 1.000000000000000e+00, };
static double ub1[] = {5.000000000000000e-01, 5.000000000000000e+00, 5.000000000000000e+00, };
static double ub2[] = {5.000000000000000e-01, 5.000000000000000e+00, 5.000000000000000e+00, };
static double ub3[] = {5.000000000000000e-01, 5.000000000000000e+00, 5.000000000000000e+00, };
static double ub4[] = {5.000000000000000e-01, 5.000000000000000e+00, 5.000000000000000e+00, };
static double ub5[] = {5.000000000000000e+00, 5.000000000000000e+00, };

static double C0[] = {};
static double C1[] = {};
static double C2[] = {};
static double C3[] = {};
static double C4[] = {};
static double C5[] = {1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, };

static double D0[] = {};
static double D1[] = {};
static double D2[] = {};
static double D3[] = {};
static double D4[] = {};
static double D5[] = {};

static double lg0[] = {};
static double lg1[] = {};
static double lg2[] = {};
static double lg3[] = {};
static double lg4[] = {};
static double lg5[] = {0.000000000000000e+00, 0.000000000000000e+00, };

static double ug0[] = {};
static double ug1[] = {};
static double ug2[] = {};
static double ug3[] = {};
static double ug4[] = {};
static double ug5[] = {0.000000000000000e+00, 0.000000000000000e+00, };

static double Zl0[] = {};
static double Zl1[] = {};
static double Zl2[] = {};
static double Zl3[] = {};
static double Zl4[] = {};
static double Zl5[] = {1.000000000000000e+06, 1.000000000000000e+06, };

static double Zu0[] = {};
static double Zu1[] = {};
static double Zu2[] = {};
static double Zu3[] = {};
static double Zu4[] = {};
static double Zu5[] = {1.000000000000000e+06, 1.000000000000000e+06, };

static double zl0[] = {};
static double zl1[] = {};
static double zl2[] = {};
static double zl3[] = {};
static double zl4[] = {};
static double zl5[] = {0.000000000000000e+00, 0.000000000000000e+00, };

static double zu0[] = {};
static double zu1[] = {};
static double zu2[] = {};
static double zu3[] = {};
static double zu4[] = {};
static double zu5[] = {0.000000000000000e+00, 0.000000000000000e+00, };

static int idxs0[] = {};
static int idxs1[] = {};
static int idxs2[] = {};
static int idxs3[] = {};
static int idxs4[] = {};
static int idxs5[] = {2, 3, };

static double lls0[] = {};
static double lls1[] = {};
static double lls2[] = {};
static double lls3[] = {};
static double lls4[] = {};
static double lls5[] = {0.000000000000000e+00, 0.000000000000000e+00, };

static double lus0[] = {};
static double lus1[] = {};
static double lus2[] = {};
static double lus3[] = {};
static double lus4[] = {};
static double lus5[] = {0.000000000000000e+00, 0.000000000000000e+00, };

static double u_guess0[] = {0.000000000000000e+00, };
static double u_guess1[] = {0.000000000000000e+00, };
static double u_guess2[] = {0.000000000000000e+00, };
static double u_guess3[] = {0.000000000000000e+00, };
static double u_guess4[] = {0.000000000000000e+00, };
static double u_guess5[] = {};

static double x_guess0[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double x_guess1[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double x_guess2[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double x_guess3[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double x_guess4[] = {0.000000000000000e+00, 0.000000000000000e+00, };
static double x_guess5[] = {0.000000000000000e+00, 0.000000000000000e+00, };

static double sl_guess0[] = {};
static double sl_guess1[] = {};
static double sl_guess2[] = {};
static double sl_guess3[] = {};
static double sl_guess4[] = {};
static double sl_guess5[] = {1.000000000000000e+00, 1.000000000000000e+00, };

static double su_guess0[] = {};
static double su_guess1[] = {};
static double su_guess2[] = {};
static double su_guess3[] = {};
static double su_guess4[] = {};
static double su_guess5[] = {1.000000000000000e+00, 1.000000000000000e+00, };

static double *AA[] = {A0, A1, A2, A3, A4};
static double *BB[] = {B0, B1, B2, B3, B4};
static double *bb[] = {b0, b1, b2, b3, b4};
static double *QQ[] = {Q0, Q1, Q2, Q3, Q4, Q5};
static double *RR[] = {R0, R1, R2, R3, R4, R5};
static double *SS[] = {S0, S1, S2, S3, S4, S5};
static double *qq[] = {q0, q1, q2, q3, q4, q5};
static double *rr[] = {r0, r1, r2, r3, r4, r5};
static int *iidxb[] = {idxb0, idxb1, idxb2, idxb3, idxb4, idxb5};
static double *llb[] = {lb0, lb1, lb2, lb3, lb4, lb5};
static double *uub[] = {ub0, ub1, ub2, ub3, ub4, ub5};
static double *CC[] = {C0, C1, C2, C3, C4, C5};
static double *DD[] = {D0, D1, D2, D3, D4, D5};
static double *llg[] = {lg0, lg1, lg2, lg3, lg4, lg5};
static double *uug[] = {ug0, ug1, ug2, ug3, ug4, ug5};
static double *ZZl[] = {Zl0, Zl1, Zl2, Zl3, Zl4, Zl5};
static double *ZZu[] = {Zu0, Zu1, Zu2, Zu3, Zu4, Zu5};
static double *zzl[] = {zl0, zl1, zl2, zl3, zl4, zl5};
static double *zzu[] = {zu0, zu1, zu2, zu3, zu4, zu5};
static int *iidxs[] = {idxs0, idxs1, idxs2, idxs3, idxs4, idxs5};
static double *llls[] = {lls0, lls1, lls2, lls3, lls4, lls5};
static double *llus[] = {lus0, lus1, lus2, lus3, lus4, lus5};
static double *uu_guess[] = {u_guess0, u_guess1, u_guess2, u_guess3, u_guess4, u_guess5};
static double *xx_guess[] = {x_guess0, x_guess1, x_guess2, x_guess3, x_guess4, x_guess5};
static double *ssl_guess[] = {sl_guess0, sl_guess1, sl_guess2, sl_guess3, sl_guess4, sl_guess5};
static double *ssu_guess[] = {su_guess0, su_guess1, su_guess2, su_guess3, su_guess4, su_guess5};

int *nu = nnu;
int *nx = nnx;
int *nbu = nnbu;
int *nbx = nnbx;
int *ng = nng;
int *ns = nns;
double **hA = AA;
double **hB = BB;
double **hb = bb;
double **hQ = QQ;
double **hR = RR;
double **hS = SS;
double **hq = qq;
double **hr = rr;
int **hidxb = iidxb;
double **hlb = llb;
double **hub = uub;
double **hC = CC;
double **hD = DD;
double **hlg = llg;
double **hug = uug;
double **hZl = ZZl;
double **hZu = ZZu;
double **hzl = zzl;
double **hzu = zzu;
int **hidxs = iidxs;
double **hlls = llls;
double **hlus = llus;
double **hu_guess = uu_guess;
double **hx_guess = xx_guess;
double **hsl_guess = ssl_guess;
double **hsu_guess = ssu_guess;
