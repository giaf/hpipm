function codegen_ocp_qp_data( dims, qp );


% filename
filename = "ocp_qp_data.c";
% open file
fid = fopen(filename, "w");


% extract dims
N = dims.N;
nx = dims.nx;
nu = dims.nu;
nbx = dims.nbx;
nbu = dims.nbu;
ng = dims.ng;
nsbx = dims.nsbx;
nsbu = dims.nsbu;
nsg = dims.nsg;


% N
fprintf(fid, "int N = %d;\n", N);
% nx
fprintf(fid, "static int nnx[] = {");
for ii=1:N
	fprintf(fid, "%d, ", nx(ii));
end
fprintf(fid, "%d};\n", nx(ii+1));
% nu
fprintf(fid, "static int nnu[] = {");
for ii=1:N
	fprintf(fid, "%d, ", nu(ii));
end
fprintf(fid, "%d};\n", nu(ii+1));
% nbx
fprintf(fid, "static int nnbx[] = {");
for ii=1:N
	fprintf(fid, "%d, ", nbx(ii));
end
fprintf(fid, "%d};\n", nbx(ii+1));
% nbu
fprintf(fid, "static int nnbu[] = {");
for ii=1:N
	fprintf(fid, "%d, ", nbu(ii));
end
fprintf(fid, "%d};\n", nbu(ii+1));
% ng
fprintf(fid, "static int nng[] = {");
for ii=1:N
	fprintf(fid, "%d, ", ng(ii));
end
fprintf(fid, "%d};\n", ng(ii+1));
% nsbx
fprintf(fid, "static int nnsbx[] = {");
for ii=1:N
	fprintf(fid, "%d, ", nsbx(ii));
end
fprintf(fid, "%d};\n", nsbx(ii+1));
% nsbu
fprintf(fid, "static int nnsbu[] = {");
for ii=1:N
	fprintf(fid, "%d, ", nsbu(ii));
end
fprintf(fid, "%d};\n", nsbu(ii+1));
% nsg
fprintf(fid, "static int nnsg[] = {");
for ii=1:N
	fprintf(fid, "%d, ", nsg(ii));
end
fprintf(fid, "%d};\n", nsg(N+1));
%
fprintf(fid, "\n");


% A
for ii=1:N
	fprintf(fid, "static double A%d[] = {", ii-1);
	jj = 0;
	jj_max = nx(ii+1)*nx(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.A{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.A{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% B
for ii=1:N
	fprintf(fid, "static double B%d[] = {", ii-1);
	jj = 0;
	jj_max = nx(ii+1)*nu(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.B{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.B{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% b
for ii=1:N
	fprintf(fid, "static double b%d[] = {", ii-1);
	jj = 0;
	jj_max = nx(ii+1);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.b{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.b{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% Q
for ii=1:N+1
	fprintf(fid, "static double Q%d[] = {", ii-1);
	jj = 0;
	jj_max = nx(ii)*nx(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.Q{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.Q{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% R
for ii=1:N+1
	fprintf(fid, "static double R%d[] = {", ii-1);
	jj = 0;
	jj_max = nu(ii)*nu(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.R{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.R{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% S
for ii=1:N+1
	fprintf(fid, "static double S%d[] = {", ii-1);
	jj = 0;
	jj_max = nu(ii)*nx(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.S{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.S{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% q
for ii=1:N+1
	fprintf(fid, "static double q%d[] = {", ii-1);
	jj = 0;
	jj_max = nx(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.q{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.q{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% r
for ii=1:N+1
	fprintf(fid, "static double r%d[] = {", ii-1);
	jj = 0;
	jj_max = nu(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.r{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.r{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% idxb
for ii=1:N+1
	fprintf(fid, "static int idxb%d[] = {", ii-1);
	jj = 0;
	jj_max = nbu(ii);
	for jj=1:jj_max-1
		kk0 = 0;
		for kk=1:nx(ii)
			if kk0==0 && qp.Ju{ii}(jj,kk)!=0
				kk0 = kk;
				fprintf(fid, "%d, ", kk-1);
			end
		end
	end
	if jj_max>0
		kk0 = 0;
		for kk=1:nx(ii)
			if kk0==0 && qp.Ju{ii}(jj+1,kk)!=0
				kk0 = kk;
				if nbx(ii)>0
					fprintf(fid, "%d, ", kk-1);
				else
					fprintf(fid, "%d", kk-1);
				end
			end
		end
	end
	jj = 0;
	jj_max = nbx(ii);
	for jj=1:jj_max-1
		kk0 = 0;
		for kk=1:nx(ii)
			if kk0==0 && qp.Jx{ii}(jj,kk)!=0
				kk0 = kk;
				fprintf(fid, "%d, ", nu(ii)+kk-1);
			end
		end
	end
	if jj_max>0
		kk0 = 0;
		for kk=1:nx(ii)
			if kk0==0 && qp.Jx{ii}(jj+1,kk)!=0
				kk0 = kk;
				fprintf(fid, "%d", nu(ii)+kk-1);
			end
		end
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% lb
for ii=1:N+1
	fprintf(fid, "static double lb%d[] = {", ii-1);
	jj = 0;
	jj_max = nbu(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.lu{ii}(:)(jj));
	end
	if jj_max>0
		if nbx(ii)>0
			fprintf(fid, "%e, ", qp.lu{ii}(:)(jj+1));
		else
			fprintf(fid, "%e", qp.lu{ii}(:)(jj+1));
		end
	end
	jj = 0;
	jj_max = nbx(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.lx{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.lx{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");
% ub
for ii=1:N+1
	fprintf(fid, "static double ub%d[] = {", ii-1);
	jj = 0;
	jj_max = nbu(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.uu{ii}(:)(jj));
	end
	if jj_max>0
		if nbx(ii)>0
			fprintf(fid, "%e, ", qp.uu{ii}(:)(jj+1));
		else
			fprintf(fid, "%e", qp.uu{ii}(:)(jj+1));
		end
	end
	jj = 0;
	jj_max = nbx(ii);
	for jj=1:jj_max-1
		fprintf(fid, "%e, ", qp.ux{ii}(:)(jj));
	end
	if jj_max>0
		fprintf(fid, "%e", qp.ux{ii}(:)(jj+1));
	end
	fprintf(fid, "};\n");
end
fprintf(fid, "\n");



% AA
fprintf(fid, "static double *AA[] = {");
for ii=1:N-1
	fprintf(fid, "A%d, ", ii-1);
end
fprintf(fid, "A%d};\n", ii);
% BB
fprintf(fid, "static double *BB[] = {");
for ii=1:N-1
	fprintf(fid, "B%d, ", ii-1);
end
fprintf(fid, "B%d};\n", ii);
% bb
fprintf(fid, "static double *bb[] = {");
for ii=1:N-1
	fprintf(fid, "b%d, ", ii-1);
end
fprintf(fid, "b%d};\n", ii);
% QQ
fprintf(fid, "static double *QQ[] = {");
for ii=1:N
	fprintf(fid, "Q%d, ", ii-1);
end
fprintf(fid, "Q%d};\n", ii);
% RR
fprintf(fid, "static double *RR[] = {");
for ii=1:N
	fprintf(fid, "R%d, ", ii-1);
end
fprintf(fid, "R%d};\n", ii);
% SS
fprintf(fid, "static double *SS[] = {");
for ii=1:N
	fprintf(fid, "S%d, ", ii-1);
end
fprintf(fid, "S%d};\n", ii);
% qq
fprintf(fid, "static double *qq[] = {");
for ii=1:N
	fprintf(fid, "q%d, ", ii-1);
end
fprintf(fid, "q%d};\n", ii);
% rr
fprintf(fid, "static double *rr[] = {");
for ii=1:N
	fprintf(fid, "r%d, ", ii-1);
end
fprintf(fid, "r%d};\n", ii);
% iidxb
fprintf(fid, "static int *iidxb[] = {");
for ii=1:N
	fprintf(fid, "idxb%d, ", ii-1);
end
fprintf(fid, "idxb%d};\n", ii);

fprintf(fid, "\n");
% llb
fprintf(fid, "static double *llb[] = {");
for ii=1:N
	fprintf(fid, "lb%d, ", ii-1);
end
fprintf(fid, "lb%d};\n", ii);

fprintf(fid, "\n");
% uub
fprintf(fid, "static double *uub[] = {");
for ii=1:N
	fprintf(fid, "ub%d, ", ii-1);
end
fprintf(fid, "ub%d};\n", ii);

fprintf(fid, "\n");



fprintf(fid, "int *nu = nnu;\n");
fprintf(fid, "int *nx = nnx;\n");
fprintf(fid, "int *nbu = nnbu;\n");
fprintf(fid, "int *nbx = nnbx;\n");
fprintf(fid, "int *ng = nng;\n");
fprintf(fid, "int *nsbu = nnsbu;\n");
fprintf(fid, "int *nsbx = nnsbx;\n");
fprintf(fid, "int *nsg = nnsg;\n");

fprintf(fid, "double **hA = AA;\n");
fprintf(fid, "double **hB = BB;\n");
fprintf(fid, "double **hb = bb;\n");
fprintf(fid, "double **hQ = QQ;\n");
fprintf(fid, "double **hR = RR;\n");
fprintf(fid, "double **hS = SS;\n");
fprintf(fid, "double **hq = qq;\n");
fprintf(fid, "double **hr = rr;\n");
fprintf(fid, "int **hidxb = iidxb;\n");
fprintf(fid, "double **hlb = llb;\n");
fprintf(fid, "double **hub = uub;\n");
fprintf(fid, "double **hC;// = CC;\n");
fprintf(fid, "double **hD;// = DD;\n");
fprintf(fid, "double **hlg;// = llg;\n");
fprintf(fid, "double **hug;// = uug;\n");
fprintf(fid, "double **hZl;// = ZZl;\n");
fprintf(fid, "double **hZu;// = ZZu;\n");
fprintf(fid, "double **hzl;// = zzl;\n");
fprintf(fid, "double **hzu;// = zzu;\n");
fprintf(fid, "int **hidxs;// = iidxs;\n");
fprintf(fid, "double **hlls;// = llls;\n");
fprintf(fid, "double **huls;// = uuls;\n");



% close file
fclose(fid);


return;
