/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* The 2-Clause BSD License                                                                        *
*                                                                                                 *
* Redistribution and use in source and binary forms, with or without                              *
* modification, are permitted provided that the following conditions are met:                     *
*                                                                                                 *
* 1. Redistributions of source code must retain the above copyright notice, this                  *
*    list of conditions and the following disclaimer.                                             *
* 2. Redistributions in binary form must reproduce the above copyright notice,                    *
*    this list of conditions and the following disclaimer in the documentation                    *
*    and/or other materials provided with the distribution.                                       *
*                                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 *
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



void DENSE_QP_DIM_PRINT(struct DENSE_QP_DIM *qp_dim)
	{
	int ii;

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int ng = qp_dim->ng;
	int nsb = qp_dim->nsb;
	int nsg = qp_dim->nsg;
	int ns = qp_dim->ns;

	printf("nv = %d\n\n", nv);
	printf("ne = %d\n\n", ne);
	printf("nb = %d\n\n", nb);
	printf("ng = %d\n\n", ng);
	printf("nsb = %d\n\n", nsb);
	printf("nsg = %d\n\n", nsg);
	printf("ns = %d\n\n", ns);

	return;
	}



void DENSE_QP_DIM_CODEGEN(char *file_name, char *mode, struct DENSE_QP_DIM *qp_dim)
	{
	int ii;

	FILE *file = fopen(file_name, mode);

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int ng = qp_dim->ng;
	int nsb = qp_dim->nsb;
	int nsg = qp_dim->nsg;
	int ns = qp_dim->ns;

	fprintf(file, "/***************\n* dim\n***************/\n");

	// nv
	fprintf(file, "/* nv */\n");
	fprintf(file, "int nv = %d;\n", nv);
	// ne
	fprintf(file, "/* ne */\n");
	fprintf(file, "int ne = %d;\n", ne);
	// nb
	fprintf(file, "/* nb */\n");
	fprintf(file, "int nb = %d;\n", nb);
	// ng
	fprintf(file, "/* ng */\n");
	fprintf(file, "int ng = %d;\n", ng);
	// nsb
	fprintf(file, "/* nsb */\n");
	fprintf(file, "int nsb = %d;\n", nsb);
	// nsg
	fprintf(file, "/* nsg */\n");
	fprintf(file, "int nsg = %d;\n", nsg);
	// ns
	fprintf(file, "/* ns */\n");
	fprintf(file, "int ns = %d;\n", ns);

	fclose(file);

	return;
	}



void DENSE_QP_PRINT(struct DENSE_QP_DIM *qp_dim, struct DENSE_QP *qp)
	{
	int ii;

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int ng = qp_dim->ng;
	int ns = qp_dim->ns;

	printf("H = \n");
	BLASFEO_PRINT_MAT(nv, nv, qp->Hv, 0, 0);

	printf("A = \n");
	BLASFEO_PRINT_MAT(ne, nv, qp->A, 0, 0);

	printf("Ct = \n");
	BLASFEO_PRINT_MAT(nv, ng, qp->Ct, 0, 0);

	printf("idxb = \n");
	int_print_mat(1, nb, qp->idxb, 1);

	printf("gz = \n");
	BLASFEO_PRINT_TRAN_VEC(nv+2*ns, qp->gz, 0);

	printf("b = \n");
	BLASFEO_PRINT_TRAN_VEC(ne, qp->b, 0);

	printf("d = \n");
	BLASFEO_PRINT_TRAN_VEC(2*nb+2*ng+2*ns, qp->d, 0);

	printf("d_mask = \n");
	BLASFEO_PRINT_TRAN_VEC(2*nb+2*ng+2*ns, qp->d_mask, 0);

	printf("m = \n");
	BLASFEO_PRINT_TRAN_VEC(2*nb+2*ng+2*ns, qp->m, 0);

	printf("Z = \n");
	BLASFEO_PRINT_TRAN_VEC(2*ns, qp->Z, 0);

	printf("idxs_rev = \n");
	int_print_mat(1, nb+ng, qp->idxs_rev, 1);

	return;
	}



void DENSE_QP_CODEGEN(char *file_name, char *mode, struct DENSE_QP_DIM *qp_dim, struct DENSE_QP *qp)
	{
	int ii, jj;

	FILE *file = fopen(file_name, mode);

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int ng = qp_dim->ng;
	int ns = qp_dim->ns;

	fprintf(file, "/***************\n* qp\n***************/\n");

	// H
	fprintf(file, "/* H */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double HH[] = {");
#else
	fprintf(file, "static float HH[] = {");
#endif
	for(jj=0; jj<nv; jj++)
		{
		for(ii=0; ii<nv; ii++)
			{
#ifdef DOUBLE_PRECISION
			fprintf(file, "%18.15e, ", BLASFEO_DMATEL(qp->Hv, ii, jj));
#else
			fprintf(file, "%18.15e, ", BLASFEO_SMATEL(qp->Hv, ii, jj));
#endif
			}
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *H = HH;\n");
#else
	fprintf(file, "float *H = HH;\n");
#endif

	// A
	fprintf(file, "/* A */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double AA[] = {");
#else
	fprintf(file, "static float AA[] = {");
#endif
	for(jj=0; jj<nv; jj++)
		{
		for(ii=0; ii<ne; ii++)
			{
#ifdef DOUBLE_PRECISION
			fprintf(file, "%18.15e, ", BLASFEO_DMATEL(qp->A, ii, jj));
#else
			fprintf(file, "%18.15e, ", BLASFEO_SMATEL(qp->A, ii, jj));
#endif
			}
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *A = AA;\n");
#else
	fprintf(file, "float *A = AA;\n");
#endif

	// C
	fprintf(file, "/* C */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double CC[] = {");
#else
	fprintf(file, "static float CC[] = {");
#endif
	for(ii=0; ii<nv; ii++)
		{
		for(jj=0; jj<ng; jj++)
			{
#ifdef DOUBLE_PRECISION
			fprintf(file, "%18.15e, ", BLASFEO_DMATEL(qp->Ct, ii, jj));
#else
			fprintf(file, "%18.15e, ", BLASFEO_SMATEL(qp->Ct, ii, jj));
#endif
			}
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *C = CC;\n");
#else
	fprintf(file, "float *C = CC;\n");
#endif

	// idxb
	fprintf(file, "/* idxb */\n");
	fprintf(file, "static int iidxb[] = {");
	for(ii=0; ii<nb; ii++)
		fprintf(file, "%d, ", qp->idxb[ii]);
	fprintf(file, "};\n");
	fprintf(file, "int *idxb = iidxb;\n");

	// g
	fprintf(file, "/* g */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double gg[] = {");
#else
	fprintf(file, "static float gg[] = {");
#endif
	for(ii=0; ii<nv; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->gz, ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->gz, ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *g = gg;\n");
#else
	fprintf(file, "float *g = gg;\n");
#endif

	// zl
	fprintf(file, "/* zl */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double zzl[] = {");
#else
	fprintf(file, "static float zzl[] = {");
#endif
	for(ii=0; ii<ns; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->gz, nv+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->gz, nv+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *zl = zzl;\n");
#else
	fprintf(file, "float *zl = zzl;\n");
#endif

	// zu
	fprintf(file, "/* zu */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double zzu[] = {");
#else
	fprintf(file, "static float zzu[] = {");
#endif
	for(ii=0; ii<ns; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->gz, nv+ns+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->gz, nv+ns+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *zu = zzu;\n");
#else
	fprintf(file, "float *zu = zzu;\n");
#endif

	// b
	fprintf(file, "/* b */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double bb[] = {");
#else
	fprintf(file, "static float bb[] = {");
#endif
	for(ii=0; ii<ne; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->b, ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->b, ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *b = bb;\n");
#else
	fprintf(file, "float *b = bb;\n");
#endif

	// lb
	fprintf(file, "/* lb */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double llb[] = {");
#else
	fprintf(file, "static float llb[] = {");
#endif
	for(ii=0; ii<nb; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d, ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d, ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *lb = llb;\n");
#else
	fprintf(file, "float *lb = llb;\n");
#endif

	// lb_mask
	fprintf(file, "/* lb_mask */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double llb_mask[] = {");
#else
	fprintf(file, "static float llb_mask[] = {");
#endif
	for(ii=0; ii<nb; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d_mask, ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d_mask, ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *lb_mask = llb_mask;\n");
#else
	fprintf(file, "float *lb_mask = llb_mask;\n");
#endif

	// ub
	fprintf(file, "/* ub */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double uub[] = {");
#else
	fprintf(file, "static float uub[] = {");
#endif
	for(ii=0; ii<nb; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", -BLASFEO_DVECEL(qp->d, nb+ng+ii));
#else
		fprintf(file, "%18.15e, ", -BLASFEO_SVECEL(qp->d, nb+ng+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *ub = uub;\n");
#else
	fprintf(file, "float *ub = uub;\n");
#endif

	// ub_mask
	fprintf(file, "/* ub_mask */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double uub_mask[] = {");
#else
	fprintf(file, "static float uub_mask[] = {");
#endif
	for(ii=0; ii<nb; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d_mask, nb+ng+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d_mask, nb+ng+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *ub_mask = uub_mask;\n");
#else
	fprintf(file, "float *ub_mask = uub_mask;\n");
#endif

	// lg
	fprintf(file, "/* lg */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double llg[] = {");
#else
	fprintf(file, "static float llg[] = {");
#endif
	for(ii=0; ii<ng; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d, nb+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d, nb+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *lg = llg;\n");
#else
	fprintf(file, "float *lg = llg;\n");
#endif

	// lg_mask
	fprintf(file, "/* lg_mask */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double llg_mask[] = {");
#else
	fprintf(file, "static float llg_mask[] = {");
#endif
	for(ii=0; ii<ng; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d_mask, nb+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d_mask, nb+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *lg_mask = llg_mask;\n");
#else
	fprintf(file, "float *lg_mask = llg_mask;\n");
#endif

	// ug
	fprintf(file, "/* ug */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double uug[] = {");
#else
	fprintf(file, "static float uug[] = {");
#endif
	for(ii=0; ii<ng; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", -BLASFEO_DVECEL(qp->d, 2*nb+ng+ii));
#else
		fprintf(file, "%18.15e, ", -BLASFEO_SVECEL(qp->d, 2*nb+ng+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *ug = uug;\n");
#else
	fprintf(file, "float *ug = uug;\n");
#endif

	// ug_mask
	fprintf(file, "/* ug_mask */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double uug_mask[] = {");
#else
	fprintf(file, "static float uug_mask[] = {");
#endif
	for(ii=0; ii<ng; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d_mask, 2*nb+ng+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d_mask, 2*nb+ng+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *ug_mask = uug_mask;\n");
#else
	fprintf(file, "float *ug_mask = uug_mask;\n");
#endif

	// lls
	fprintf(file, "/* lls */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double llls[] = {");
#else
	fprintf(file, "static float llls[] = {");
#endif
	for(ii=0; ii<ns; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d, 2*nb+2*ng+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d, 2*nb+2*ng+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *lls = llls;\n");
#else
	fprintf(file, "float *lls = llls;\n");
#endif

	// lls_mask
	fprintf(file, "/* lls_mask */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double llls_mask[] = {");
#else
	fprintf(file, "static float llls_mask[] = {");
#endif
	for(ii=0; ii<ns; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d_mask, 2*nb+2*ng+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d_mask, 2*nb+2*ng+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *lls_mask = llls_mask;\n");
#else
	fprintf(file, "float *lls_mask = llls_mask;\n");
#endif

	// lus
	fprintf(file, "/* lus */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double llus[] = {");
#else
	fprintf(file, "static float llus[] = {");
#endif
	for(ii=0; ii<ns; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d, 2*nb+2*ng+ns+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d, 2*nb+2*ng+ns+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *lus = llus;\n");
#else
	fprintf(file, "float *lus = llus;\n");
#endif

	// lus_mask
	fprintf(file, "/* lus_mask */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double llus_mask[] = {");
#else
	fprintf(file, "static float llus_mask[] = {");
#endif
	for(ii=0; ii<ns; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->d_mask, 2*nb+2*ng+ns+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->d_mask, 2*nb+2*ng+ns+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *lus_mask = llus_mask;\n");
#else
	fprintf(file, "float *lus_mask = llus_mask;\n");
#endif

	printf("Z = \n");
	BLASFEO_PRINT_TRAN_VEC(2*ns, qp->Z, 0);
	// zl
	fprintf(file, "/* Zl */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double ZZl[] = {");
#else
	fprintf(file, "static float ZZl[] = {");
#endif
	for(ii=0; ii<ns; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->Z, ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->Z, ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *Zl = ZZl;\n");
#else
	fprintf(file, "float *Zl = ZZl;\n");
#endif

	// Zu
	fprintf(file, "/* Zu */\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "static double ZZu[] = {");
#else
	fprintf(file, "static float ZZu[] = {");
#endif
	for(ii=0; ii<ns; ii++)
		{
#ifdef DOUBLE_PRECISION
		fprintf(file, "%18.15e, ", BLASFEO_DVECEL(qp->Z, ns+ii));
#else
		fprintf(file, "%18.15e, ", BLASFEO_SVECEL(qp->Z, ns+ii));
#endif
		}
	fprintf(file, "};\n");
#ifdef DOUBLE_PRECISION
	fprintf(file, "double *Zu = ZZu;\n");
#else
	fprintf(file, "float *Zu = ZZu;\n");
#endif

	// idxs_rev
	fprintf(file, "/* idxs_rev */\n");
	fprintf(file, "static int iidxs_rev[] = {");
	for(ii=0; ii<nb+ng; ii++)
		fprintf(file, "%d, ", qp->idxs_rev[ii]);
	fprintf(file, "};\n");
	fprintf(file, "int *idxs_rev = iidxs_rev;\n");

	fclose(file);

	return;
	}



void DENSE_QP_SOL_PRINT(struct DENSE_QP_DIM *qp_dim, struct DENSE_QP_SOL *qp_sol)
	{
	int ii;

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int ng = qp_dim->ng;
	int ns = qp_dim->ns;

	printf("v = \n");
	BLASFEO_PRINT_TRAN_VEC(nv+2*ns, qp_sol->v, 0);

	printf("pi = \n");
	BLASFEO_PRINT_TRAN_VEC(ne, qp_sol->pi, 0);

	printf("lam = \n");
	BLASFEO_PRINT_TRAN_VEC(2*nb+2*ng+2*ns, qp_sol->lam, 0);

	printf("t = \n");
	BLASFEO_PRINT_TRAN_VEC(2*nb+2*ng+2*ns, qp_sol->t, 0);

	return;
	}



void DENSE_QP_RES_PRINT(struct DENSE_QP_DIM *qp_dim, struct DENSE_QP_RES *qp_res)
	{
	int ii;

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int ng = qp_dim->ng;
	int ns = qp_dim->ns;

	printf("res_g = \n");
	BLASFEO_PRINT_TRAN_VEC(nv+2*ns, qp_res->res_g, 0);

	printf("res_b = \n");
	BLASFEO_PRINT_TRAN_VEC(ne, qp_res->res_b, 0);

	printf("res_d = \n");
	BLASFEO_PRINT_TRAN_VEC(2*nb+2*ng+2*ns, qp_res->res_d, 0);

	printf("res_m = \n");
	BLASFEO_PRINT_TRAN_VEC(2*nb+2*ng+2*ns, qp_res->res_m, 0);

	return;
	}



void DENSE_QP_IPM_ARG_PRINT(struct DENSE_QP_DIM *qp_dim, struct DENSE_QP_IPM_ARG *arg)
	{
	int ii;

	// mode
	printf("/* mode */\n");
	printf("int mode = %d;\n", arg->mode);
	// iter_max
	printf("/* iter_max */\n");
	printf("int iter_max = %d;\n", arg->iter_max);
	// alpha_min
	printf("/* alpha_min */\n");
	printf("double alpha_min = %18.15e;\n", arg->alpha_min);
	// mu0
	printf("/* mu0 */\n");
	printf("double mu0 = %18.15e;\n", arg->mu0);
	// tol_stat
	printf("/* tol_stat */\n");
	printf("double tol_stat = %18.15e;\n", arg->res_g_max);
	// tol_eq
	printf("/* tol_eq */\n");
	printf("double tol_eq = %18.15e;\n", arg->res_b_max);
	// tol_ineq
	printf("/* tol_ineq */\n");
	printf("double tol_ineq = %18.15e;\n", arg->res_d_max);
	// tol_comp
	printf("/* tol_comp */\n");
	printf("double tol_comp = %18.15e;\n", arg->res_m_max);
	// reg_prim
	printf("/* reg_prim */\n");
	printf("double reg_prim = %18.15e;\n", arg->reg_prim);
	// reg_dual
	printf("/* reg_dual */\n");
	printf("double reg_dual = %18.15e;\n", arg->reg_dual);
	// warm_start
	printf("/* warm_start */\n");
	printf("int warm_start = %d;\n", arg->warm_start);
	// pred_corr
	printf("/* pred_corr */\n");
	printf("int pred_corr = %d;\n", arg->pred_corr);
	// split_step
	printf("/* split_step */\n");
	printf("int split_step = %d;\n", arg->split_step);

	return;
	}



void DENSE_QP_IPM_ARG_CODEGEN(char *file_name, char *mode, struct DENSE_QP_DIM *qp_dim, struct DENSE_QP_IPM_ARG *arg)
	{
	int ii;

	FILE *file = fopen(file_name, mode);

	fprintf(file, "/***************\n* arg\n***************/\n");

	// mode
	fprintf(file, "/* mode */\n");
	fprintf(file, "int mode = %d;\n", arg->mode);
	// iter_max
	fprintf(file, "/* iter_max */\n");
	fprintf(file, "int iter_max = %d;\n", arg->iter_max);
	// alpha_min
	fprintf(file, "/* alpha_min */\n");
	fprintf(file, "double alpha_min = %18.15e;\n", arg->alpha_min);
	// mu0
	fprintf(file, "/* mu0 */\n");
	fprintf(file, "double mu0 = %18.15e;\n", arg->mu0);
	// tol_stat
	fprintf(file, "/* tol_stat */\n");
	fprintf(file, "double tol_stat = %18.15e;\n", arg->res_g_max);
	// tol_eq
	fprintf(file, "/* tol_eq */\n");
	fprintf(file, "double tol_eq = %18.15e;\n", arg->res_b_max);
	// tol_ineq
	fprintf(file, "/* tol_ineq */\n");
	fprintf(file, "double tol_ineq = %18.15e;\n", arg->res_d_max);
	// tol_comp
	fprintf(file, "/* tol_comp */\n");
	fprintf(file, "double tol_comp = %18.15e;\n", arg->res_m_max);
	// reg_prim
	fprintf(file, "/* reg_prim */\n");
	fprintf(file, "double reg_prim = %18.15e;\n", arg->reg_prim);
	// reg_dual
	fprintf(file, "/* reg_dual */\n");
	fprintf(file, "double reg_dual = %18.15e;\n", arg->reg_dual);
	// warm_start
	fprintf(file, "/* warm_start */\n");
	fprintf(file, "int warm_start = %d;\n", arg->warm_start);
	// pred_corr
	fprintf(file, "/* pred_corr */\n");
	fprintf(file, "int pred_corr = %d;\n", arg->pred_corr);
	// split_step
	fprintf(file, "/* split_step */\n");
	fprintf(file, "int split_step = %d;\n", arg->split_step);

	fclose(file);

	return;
	}




