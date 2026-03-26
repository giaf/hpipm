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



void COMPUTE_GAMMA_GAMMA_QP(REAL *res_d, REAL *res_m, struct CORE_QP_IPM_WORKSPACE *cws)
	{

	int nc = cws->nc;

	REAL *lam = cws->lam;
	REAL *t = cws->t;
//	REAL *res_d = cws->res_d; // TODO rename d ???
//	REAL *res_m = cws->res_m; // TODO rename m ???
	REAL *t_inv = cws->t_inv;
	REAL *Gamma = cws->Gamma;
	REAL *gamma = cws->gamma;
	REAL lam_min = cws->lam_min;
	REAL t_min = cws->t_min;
	REAL t_min_inv = cws->t_min_inv;
	REAL lam0, t0, t_inv_tmp, lam_tmp;

	// local variables
	int ii;

//printf("\ngamma gamma\n");
	if(cws->t_lam_min==1)
		{
		for(ii=0; ii<nc; ii++)
			{
			lam0 = lam[ii];
			t0 = t[ii];
			t_inv[ii] = 1.0/t0;
			t_inv_tmp = t0<t_min ? t_min_inv : t_inv[ii];
			lam_tmp = lam0<lam_min ? lam_min : lam0;
			Gamma[ii] = t_inv_tmp*lam_tmp;
			gamma[ii] = t_inv[ii]*(res_m[ii]-lam0*res_d[ii]);
			}
		}
	else
		{
		for(ii=0; ii<nc; ii++)
			{
			lam0 = lam[ii];
			t0 = t[ii];
			t_inv[ii] = 1.0/t0;
			Gamma[ii] = t_inv[ii]*lam0;
			gamma[ii] = t_inv[ii]*(res_m[ii]-lam0*res_d[ii]);
			}
		}

	return;

	}



void COMPUTE_GGAMMA_QP(struct CORE_QP_IPM_WORKSPACE *cws)
	{

	int nc = cws->nc;

	REAL *lam = cws->lam;
	REAL *t = cws->t;
	REAL *t_inv = cws->t_inv;
	REAL *Gamma = cws->Gamma;
	REAL lam_min = cws->lam_min;
	REAL t_min = cws->t_min;
	REAL t_min_inv = cws->t_min_inv;
	REAL lam0, t0, t_inv_tmp, lam_tmp;

	// local variables
	int ii;

//printf("\nGamma\n");
	if(cws->t_lam_min==1)
		{
		for(ii=0; ii<nc; ii++)
			{
			lam0 = lam[ii];
			t0 = t[ii];
			t_inv[ii] = 1.0/t0;
			t_inv_tmp = t0<t_min ? t_min_inv : t_inv[ii];
			lam_tmp = lam0<lam_min ? lam_min : lam0;
			Gamma[ii] = t_inv_tmp*lam_tmp;
			}
		}
	else
		{
		for(ii=0; ii<nc; ii++)
			{
			lam0 = lam[ii];
			t0 = t[ii];
			t_inv[ii] = 1.0/t0;
			Gamma[ii] = t_inv[ii]*lam0;
			}
		}

	return;

	}



void COMPUTE_GAMMA_QP(REAL *res_d, REAL *res_m, struct CORE_QP_IPM_WORKSPACE *cws)
	{

	int nc = cws->nc;

	REAL *lam = cws->lam;
//	REAL *res_m = cws->res_m;
//	REAL *res_d = cws->res_d;
	REAL *t_inv = cws->t_inv;
	REAL *gamma = cws->gamma;
	REAL lam_min = cws->lam_min;
	REAL lam0;

	// local variables
	int ii;

	for(ii=0; ii<nc; ii++)
		{
		lam0 = lam[ii];
//		lam0 = lam0<lam_min ? lam_min : lam0;
		gamma[ii] = t_inv[ii]*(res_m[ii]-lam0*res_d[ii]);
		}

	return;

	}

void COMPUTE_LAM_T_QP(REAL *res_d, REAL *res_m, REAL *dlam, REAL *dt, struct CORE_QP_IPM_WORKSPACE *cws)
	{

	int nc = cws->nc;

	REAL *lam = cws->lam;
	REAL *t_inv = cws->t_inv;
	REAL lam_min = cws->lam_min;
	REAL lam0;

	// local variables
	int ii;

	for(ii=0; ii<nc; ii++)
		{
		lam0 = lam[ii];
//		lam0 = lam0<lam_min ? lam_min : lam0;
		dlam[ii] = - t_inv[ii] * (res_m[ii] + (lam0*dt[ii]) - (lam0*res_d[ii]));
		dt[ii] -= res_d[ii];
		// TODO compute lamda alone ???
//		dlam[ii] = - t_inv[ii] * (lam0*dt[ii] + res_m[ii]);
		}

	return;

	}



void COMPUTE_ALPHA_QP(struct CORE_QP_IPM_WORKSPACE *cws)
	{
	
	// extract workspace members
	int nc = cws->nc;

	REAL *lam = cws->lam;
	REAL *t = cws->t;
	REAL *dlam = cws->dlam;
	REAL *dt = cws->dt;

	REAL *m = cws->m;

	REAL alpha_prim = 1.0;
	REAL alpha_dual = 1.0;
	REAL alpha = 1.0;

	int m_zero = cws->m_zero;

	REAL m_safe = cws->m_safe; //0.3; // in [0.0,1.0]
	m_safe = m_safe>=0.0 ? m_safe : 0.0;
	m_safe = m_safe<=1.0 ? m_safe : 1.0;

	// local variables
	int ii;

	if(cws->split_step==1)
		{

		if(m_zero)
			{
			for(ii=0; ii<nc; ii++)
				{
				REAL lam1 = lam[ii]+alpha_dual*dlam[ii];
				REAL t1 = t[ii]+alpha_prim*dt[ii];
				if( lam1<0.0 )
					{
					alpha_dual = - lam[ii] / dlam[ii];
					lam1 = lam[ii]+alpha_dual*dlam[ii];
					}
				if( t1<0.0 )
					{
					alpha_prim = - t[ii] / dt[ii];
					t1 = t[ii]+alpha_prim*dt[ii];
					}
				}
			}
		else
			{

			// 1st pass
			for(ii=0; ii<nc; ii++)
				{

				REAL lam1 = lam[ii]+alpha_dual*dlam[ii];
				REAL t1 = t[ii]+alpha_prim*dt[ii];
				if( lam1<0.0 )
					{
					alpha_dual = - lam[ii] / dlam[ii];
					lam1 = lam[ii]+alpha_dual*dlam[ii];
					//printf("new alpha_dual %e\n", alpha_dual);
					}
				if( t1<0.0 )
					{
					alpha_prim = - t[ii] / dt[ii];
					t1 = t[ii]+alpha_prim*dt[ii];
					//printf("new alpha_prim %e\n", alpha_prim);
					}

				REAL m1 = m_safe*m[ii];

				if(lam1*t1-m1<-1e-12)
					{
					if(dlam[ii]<0)
						{
						if(dt[ii]<0)
							{
							// ignore for this pass
							//printf("\nignore for now\n");
							}
						else
							{
							alpha_dual = (m1 - lam[ii]*t1) / (dlam[ii]*t1); // XXX later thighter values of alpha_prim would change this !
							//alpha_dual = (m1 - lam[ii]*t[ii]) / (dlam[ii]*t[ii]); // be conservative and use t instead (i.e. worst case alpha_prim=0)
							}
						}
					else if(dt[ii]<0)
						{
						alpha_prim = (m1 - t[ii]*lam1) / (dt[ii]*lam1); // XXX later thighter values of alpha_dual would change this !
						//alpha_prim = (m1 - t[ii]*lam[ii]) / (dt[ii]*lam[ii]); // XXX be conservative and use lam instead (i.e. worst case alpha_dual=0)
						}
					//printf("check lam1 %e %e t1 %e %e prod %e %e\n", lam1, lam[ii]+alpha_dual*dlam[ii], t1, t[ii]+alpha_prim*dt[ii], lam1*t1-m1, (lam[ii]+alpha_dual*dlam[ii])*(t[ii]+alpha_prim*dt[ii])-m1);
					}

				}

			//printf("alpha prim %e dual %e\n", alpha_prim, alpha_dual);

			// 2nd pass, quadratic
			for(ii=0; ii<nc; ii++)
				{
				if(dlam[ii]<0 & dt[ii]<0)
					{
					REAL lam1 = lam[ii]+alpha_dual*dlam[ii];
					REAL t1 = t[ii]+alpha_prim*dt[ii];
					REAL m1 = m_safe*m[ii];
					//printf("quad prod %e\n", lam1*t1-m1);

					if(lam1*t1-m1<-1e-12)
						{
						//printf("violated prod %e\n", lam1*t1-m1);
						REAL c = lam[ii]*t[ii] - m1; // >0
						if(c>0.0)
							{
							REAL adlam = alpha_dual*dlam[ii]; // scaled with current alpha_dual
							REAL adt = alpha_prim*dt[ii]; // scaled with current alpha_prim
							REAL a = adlam*adt; // >0
							REAL b = adlam*t[ii]+lam[ii]*adt; // <0
							REAL d = b*b - 4.0*a*c; // <b*b
							//REAL r0 = 1.0;
							REAL r1 = 1.0;
							// there must be a solution between 0 (i.e. lam*t>m1) and 1 (i.e. lam1*t1<m1)
							REAL sd = sqrt(d); // <b
							REAL tmp = 0.5/a; // >0
							// it always holds 0 < r1 < 1 and r0 > r1
							//r0 = (- b + sd) * tmp; // >0
							r1 = (- b - sd) * tmp; // 0<r1<r0
							r1 = r1>0.0 ? r1 : 0.0; // if lam*t<m from previous iteration, c<0 and r1<0, so setting r1=0 triggers min_step
							//printf("\nr0 %e r1 %e\n", r0, r1);
							//printf("\nr1 %e\n", r1);
							//if(r1>0.0 && r1<1.0)
								{
								alpha_dual *= r1;
								alpha_prim *= r1;
								}
							//else
							//	{
							//	// printf("\nsomething wrong\n");
							//	}
							}
						else // lam*t<m because of initialization or previous step too long: trigger minimum step error
							{
							alpha_prim = 0.0;
							alpha_dual = 0.0;
							}
						}
					}
				}

			// 3rd pass, decreasing error in m prod
			for(ii=0; ii<nc; ii++)
				{
				REAL lam1 = lam[ii]+alpha_dual*dlam[ii];
				REAL t1 = t[ii]+alpha_prim*dt[ii];
				REAL m1 = m_safe*m[ii];

				if(lam1*t1-m1<-1e-12)
					{
					//printf("violated prod %e\n", lam1*t1-m1);
					if(dlam[ii]<0)
						{
						if(dt[ii]<0)
							{
							// not blocking in this pass
							}
						else
							{
							alpha_dual = (m1 - lam[ii]*t1) / (dlam[ii]*t1); // XXX later thighter values of alpha_prim would change this !
							}
						}
					else if(dt[ii]<0)
						{
						alpha_prim = (m1 - t[ii]*lam1) / (dt[ii]*lam1); // XXX later thighter values of alpha_dual would change this !
						}
					//printf("check lam1 %e %e t1 %e %e prod %e %e\n", lam1, lam[ii]+alpha_dual*dlam[ii], t1, t[ii]+alpha_prim*dt[ii], lam1*t1-m1, (lam[ii]+alpha_dual*dlam[ii])*(t[ii]+alpha_prim*dt[ii])-m1);
					}
				}
			}

		//printf("alpha prim %e dual %e\n\n", alpha_prim, alpha_dual);

		}
	else
		{

		if(m_zero)
			{
			for(ii=0; ii<nc; ii++)
				{
				REAL lam1 = lam[ii]+alpha*dlam[ii];
				REAL t1 = t[ii]+alpha*dt[ii];
				if( lam1<0.0 )
					{
					alpha = - lam[ii] / dlam[ii];
					lam1 = lam[ii]+alpha*dlam[ii];
					}
				if( t1<0.0 )
					{
					alpha = - t[ii] / dt[ii];
					t1 = t[ii]+alpha*dt[ii];
					}
				}
			}
		else
			{
			for(ii=0; ii<nc; ii++)
				{

				REAL lam1 = lam[ii]+alpha*dlam[ii];
				REAL t1 = t[ii]+alpha*dt[ii];
				if( lam1<0.0 )
					{
					alpha = - lam[ii] / dlam[ii];
					lam1 = lam[ii]+alpha*dlam[ii];
					//printf("new alpha (from dual) %e\n", alpha);
					}
				if( t1<0.0 )
					{
					alpha = - t[ii] / dt[ii];
					t1 = t[ii]+alpha*dt[ii];
					//printf("new alpha (from primal) %e\n", alpha);
					}

				REAL m1 = m_safe*m[ii];

				if(lam1*t1-m1<-1e-12)
					{
					// (dlam<0 & dt<0) | (dlam<0 | dt<0)
					REAL c = lam[ii]*t[ii] - m1; // >0 | >0
					if(c>0.0)
						{
						REAL a = dlam[ii]*dt[ii]; // >0 | <0
						REAL b = dlam[ii]*t[ii]+lam[ii]*dt[ii]; // <0 | ?
						REAL d = b*b - 4.0*a*c; // <b*b | >b*b
						REAL sd = sqrt(d); // <b | >b
						REAL tmp = 0.5/a; // >0 | <0
						//REAL r0 = (- b + sd) * tmp; // >0 | <0
						REAL r1 = (- b - sd) * tmp; // 0<r1<r0 | >0
						//printf("nn r0 %e r1 %e\n", r0, r1);
						//if(r1<alpha) // always true ???
						//	{
							alpha = r1;
						//	}
						//else
						//	{
						//	printf("\nimpossible\n");
						//	}
						}
					else // lam*t<m because of initialization or previous step too long: trigger minimum step error
						{
						alpha = 0.0;
						}
					}
				}
			}

		}

	if(cws->split_step==1)
		{
		cws->alpha_prim = alpha_prim;
		cws->alpha_dual = alpha_dual;
		}
	else
		{
		//REAL alpha = alpha_prim<alpha_dual ? alpha_prim : alpha_dual;
		cws->alpha_prim = alpha;
		cws->alpha_dual = alpha;
		}

	// store alpha
	//cws->alpha = alpha;

	return;

	}
	


void UPDATE_VAR_QP(struct CORE_QP_IPM_WORKSPACE *cws)
	{
	
	// extract workspace members
	int nv = cws->nv;
	int ne = cws->ne;
	int nc = cws->nc;

	REAL *v = cws->v;
	REAL *pi = cws->pi;
	REAL *lam = cws->lam;
	REAL *t = cws->t;
	REAL *v_bkp = cws->v_bkp;
	REAL *pi_bkp = cws->pi_bkp;
	REAL *lam_bkp = cws->lam_bkp;
	REAL *t_bkp = cws->t_bkp;
	REAL *dv = cws->dv;
	REAL *dpi = cws->dpi;
	REAL *dlam = cws->dlam;
	REAL *dt = cws->dt;
	//REAL alpha = cws->alpha;
	REAL alpha_prim = cws->alpha_prim;
	REAL alpha_dual = cws->alpha_dual;
	REAL lam_min = cws->lam_min;
	REAL t_min = cws->t_min;

	REAL tmp_alpha_prim, tmp_alpha_dual;
	
#if 0
	if(alpha<1.0)
		alpha *= 0.995;
#else
	REAL alpha = alpha_prim<alpha_dual ? alpha_prim : alpha_dual;
	if(alpha<1.0)
		{
		alpha_prim = alpha_prim * ((1.0-alpha_prim)*0.99 + alpha_prim*0.9999999);
		alpha_dual = alpha_dual * ((1.0-alpha_dual)*0.99 + alpha_dual*0.9999999);
		//alpha = alpha * ((1.0-alpha)*0.99 + alpha*0.9999999);
		}
#endif

	//alpha *= 0.999;
	//alpha_prim *= 0.999;
	//alpha_dual *= 0.999;

	// local variables
	int ii;

	//if(cws->split_step==0)
	//	{
	//	tmp_alpha_prim = alpha;
	//	tmp_alpha_dual = alpha;
	//	}
	//else
	//	{
		tmp_alpha_prim = alpha_prim;
		tmp_alpha_dual = alpha_dual;
	//	}

	// update v
	for(ii=0; ii<nv; ii++)
		{
		v_bkp[ii] = v[ii];
		v[ii] += tmp_alpha_prim * dv[ii];
		}

	// update pi
	for(ii=0; ii<ne; ii++)
		{
		pi_bkp[ii] = pi[ii];
		pi[ii] += tmp_alpha_dual * dpi[ii];
		}

	if(cws->t_lam_min==2)
		{
		// update lam
		for(ii=0; ii<nc; ii++)
			{
			lam_bkp[ii] = lam[ii];
			lam[ii] += tmp_alpha_dual * dlam[ii];
			lam[ii] = lam[ii]<=lam_min ? lam_min : lam[ii];
			}

		// update t
		for(ii=0; ii<nc; ii++)
			{
			t_bkp[ii] = t[ii];
			t[ii] += tmp_alpha_prim * dt[ii];
			t[ii] = t[ii]<=t_min ? t_min : t[ii];
			}
		}
	else
		{
		// update lam
		for(ii=0; ii<nc; ii++)
			{
			lam_bkp[ii] = lam[ii];
			lam[ii] += tmp_alpha_dual * dlam[ii];
			}

		// update t
		for(ii=0; ii<nc; ii++)
			{
			t_bkp[ii] = t[ii];
			t[ii] += tmp_alpha_prim * dt[ii];
			}
		}

	return;

	}



void BACKUP_VAR_QP(struct CORE_QP_IPM_WORKSPACE *cws)
	{
	
	// extract workspace members
	int nv = cws->nv;
	int ne = cws->ne;
	int nc = cws->nc;

	REAL *v = cws->v;
	REAL *pi = cws->pi;
	REAL *lam = cws->lam;
	REAL *t = cws->t;
	REAL *v_bkp = cws->v_bkp;
	REAL *pi_bkp = cws->pi_bkp;
	REAL *lam_bkp = cws->lam_bkp;
	REAL *t_bkp = cws->t_bkp;

	// local variables
	int ii;

	// backup v
	for(ii=0; ii<nv; ii++)
		{
		v_bkp[ii] = v[ii];
		}

	// backup pi
	for(ii=0; ii<ne; ii++)
		{
		pi_bkp[ii] = pi[ii];
		}

	// backup lam
	for(ii=0; ii<nc; ii++)
		{
		lam_bkp[ii] = lam[ii];
		}

	// backup t
	for(ii=0; ii<nc; ii++)
		{
		t_bkp[ii] = t[ii];
		}

	return;

	}



void COMPUTE_MU_AFF_QP(struct CORE_QP_IPM_WORKSPACE *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	REAL *ptr_m = cws->m;
	REAL *ptr_lam = cws->lam;
	REAL *ptr_t = cws->t;
	REAL *ptr_dlam = cws->dlam;
	REAL *ptr_dt = cws->dt;
	//REAL alpha = cws->alpha;
	REAL alpha_dual = cws->alpha_prim;
	REAL alpha_prim = cws->alpha_prim;
	// this affects the minimum value of signa !!!
//		alpha *= 0.99;

	REAL mu = 0;

	for(ii=0; ii<nc; ii++)
		{
		//mu += (ptr_lam[ii+0] + alpha*ptr_dlam[ii+0]) * (ptr_t[ii+0] + alpha*ptr_dt[ii+0]);
		//mu += fabs(- ptr_m[ii+0] + (ptr_lam[ii+0] + alpha*ptr_dlam[ii+0]) * (ptr_t[ii+0] + alpha*ptr_dt[ii+0]));
		mu += fabs(- ptr_m[ii+0] + (ptr_lam[ii+0] + alpha_dual*ptr_dlam[ii+0]) * (ptr_t[ii+0] + alpha_prim*ptr_dt[ii+0]));
		}
	
	cws->mu_aff = mu*cws->nc_mask_inv;

	return;

	}



void BACKUP_RES_M(struct CORE_QP_IPM_WORKSPACE *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	REAL *ptr_res_m = cws->res_m;
	REAL *ptr_res_m_bkp = cws->res_m_bkp;

	for(ii=0; ii<nc; ii++)
		{
		ptr_res_m_bkp[ii+0] = ptr_res_m[ii+0];
		}

	return;

	}



// TODO sigma_mu*d_mask
void COMPUTE_CENTERING_CORRECTION_QP(struct CORE_QP_IPM_WORKSPACE *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	REAL *ptr_dlam = cws->dlam;
	REAL *ptr_dt = cws->dt;
	REAL *ptr_res_m = cws->res_m;
	REAL *ptr_res_m_bkp = cws->res_m_bkp;
	REAL *weight = cws->weight;

	REAL sigma_mu = cws->sigma*cws->mu;
	sigma_mu = sigma_mu>cws->tau_min ? sigma_mu : cws->tau_min;
	cws->tau_iter = sigma_mu;

	//printf("dt*dlam:");
	//for(ii=0; ii<nc; ii++)
	//	printf(" %e", ptr_dt[ii]*ptr_dlam[ii]);
	//printf("\n");

	if(cws->use_weight)
		{
		for(ii=0; ii<nc; ii++)
			{
			ptr_res_m[ii+0] = ptr_res_m_bkp[ii+0] + ptr_dt[ii+0] * ptr_dlam[ii+0] - weight[ii+0]*sigma_mu;
			}
		}
	else
		{
		for(ii=0; ii<nc; ii++)
			{
			ptr_res_m[ii+0] = ptr_res_m_bkp[ii+0] + ptr_dt[ii+0] * ptr_dlam[ii+0] - sigma_mu;
			}
		}

	return;

	}



// TODO sigma_mu*d_mask
void COMPUTE_CENTERING_QP(struct CORE_QP_IPM_WORKSPACE *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	REAL *ptr_res_m = cws->res_m;
	REAL *ptr_res_m_bkp = cws->res_m_bkp;
	REAL *weight = cws->weight;
	REAL *m = cws->m;

	REAL sigma_mu = cws->sigma*cws->mu;
	sigma_mu = sigma_mu>cws->tau_min ? sigma_mu : cws->tau_min;
	cws->tau_iter = sigma_mu;

	if(cws->use_weight)
		{
		for(ii=0; ii<nc; ii++)
			{
			ptr_res_m[ii+0] = ptr_res_m_bkp[ii+0] - 0*m[ii+0] - weight[ii+0]*sigma_mu;
			}
		}
	else
		{
		for(ii=0; ii<nc; ii++)
			{
			ptr_res_m[ii+0] = ptr_res_m_bkp[ii+0] - 0*m[ii+0] - sigma_mu;
			}
		}

	return;

	}



// TODO sigma_mu*d_mask
void COMPUTE_TAU_MIN_QP(struct CORE_QP_IPM_WORKSPACE *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	REAL *ptr_res_m = cws->res_m;
	REAL *ptr_res_m_bkp = cws->res_m_bkp;
	REAL *weight = cws->weight;

	REAL *ptr_m = cws->m;

	REAL tau_min = cws->tau_min;

	if(cws->use_weight)
		{
		for(ii=0; ii<nc; ii++)
			{
			ptr_res_m[ii+0] = ptr_res_m_bkp[ii+0] - weight[ii+0]*tau_min;
			}
		}
	else
		{
		for(ii=0; ii<nc; ii++)
			{
			ptr_res_m[ii+0] = ptr_res_m_bkp[ii+0] - tau_min;
			}
		}

	return;

	}



