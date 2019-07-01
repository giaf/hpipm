clear all
close all
clc

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	disp('ERROR: env.sh has not been sourced! Before executing this example, run:');
	disp('source env.sh');
	return;
end



% compile mex interface
compile_mex_ocp_qp();


% dim
N = 5;

tic
dim = hpipm_ocp_qp_dim(N);
tmp_time = toc;
fprintf('create dim time %e\n', tmp_time);

dim.C_dim

tic
dim.set('nx', 2, 0, N);
tmp_time = toc;
fprintf('set nx time %e\n', tmp_time);
dim.set('nu', 1, 0, N-1);
dim.set('nbx', 2, 0);
dim.set('nbx', 2, 5);

% print to shell
dim.print_C_struct();
% codegen
dim.codegen('qp_data.c', 'w');


% data
A = [1, 1; 0, 1];
B = [0; 1];
%b = [0; 0]

Q = [1, 1; 0, 1];
S = [0, 0];
R = [1];
q = [1; 1];
%r = [0];

Jx = [1, 0; 0, 1];
x0 = [1; 1];


% qp
tic
qp = hpipm_ocp_qp(dim);
tmp_time = toc;
fprintf('create qp time %e\n', tmp_time);

qp.C_qp

tic
qp.set('A', A, 0, N-1);
tmp_time = toc;
fprintf('create set A time %e\n', tmp_time);
qp.set('B', B, 0, N-1);
qp.set('Q', Q, 0, N);
qp.set('S', S, 0, N-1);
qp.set('R', R, 0, N-1);
qp.set('q', q, 0, N);
%qp.set('r', r, 0, N-1);
qp.set('Jx', Jx, 0);
qp.set('lx', x0, 0);
qp.set('ux', x0, 0);
qp.set('Jx', Jx, N);

% print to shell
qp.print_C_struct();
% codegen
qp.codegen('qp_data.c', 'a');




if is_octave()
	% directly call destructor for octave 4.2.2 (ubuntu 16.04) + others ???
	if strcmp(version(), '4.2.2')
		delete(dim);
		delete(qp);
	end
end



return




% data
A = [1, 1; 0, 1];
B = [0; 1];
%b = [0; 0]

Q = [1, 1; 0, 1];
S = [0, 0];
R = [1];
q = [1; 1];
%r = [0];

Jx = [1, 0; 0, 1];
x0 = [1; 1];
Jsx = [1, 0; 0, 1];
Zl = [1e5, 0; 0, 1e5];
Zu = [1e5, 0; 0, 1e5];
zl = [1e5; 1e5];
zu = [1e5; 1e5];



% qp
tic
qp = hpipm_ocp_qp(dims);
tmp_time = toc;
fprintf('create qp time %e\n', tmp_time);

tic
qp.set('A', {A, A, A, A, A});
tmp_time = toc;
fprintf('create set A time %e\n', tmp_time);
qp.set('B', {B, B, B, B, B});
%qp.set_b({b, b, b, b, b});

qp.set('Q', {Q, Q, Q, Q, Q, Q});
qp.set('S', {S, S, S, S, S});
qp.set('R', {R, R, R, R, R});
qp.set('q', {q, q, q, q, q, q});
qp.set('q', q, 0);
%qp.set_r({r, r, r, r, r});
qp.set('Jx', Jx, 0);
qp.set('lx', x0, 0);
qp.set('ux', x0, 0);
qp.set('Jx', Jx, 5);

qp.print_C_struct()



% qp sol
tic
qp_sol = hpipm_ocp_qp_sol(dims);
tmp_time = toc;
fprintf('create qp_sol time %e\n', tmp_time);



% set up solver arg
tic
arg = hpipm_ocp_qp_solver_arg(dims);
tmp_time = toc;
fprintf('create solver arg time %e\n', tmp_time);

arg.set_mu0(1e4);
arg.set_iter_max(30);
arg.set_tol_stat(1e-4);
arg.set_tol_eq(1e-5);
arg.set_tol_ineq(1e-5);
arg.set_tol_comp(1e-5);
arg.set_reg_prim(1e-12);


% set up solver
tic
solver = hpipm_ocp_qp_solver(dims, arg);
tmp_time = toc;
fprintf('create solver time %e\n', tmp_time);



% solve qp
tic
return_flag = solver.solve(qp, qp_sol);
tmp_time = toc;
fprintf('solve time %e\n', tmp_time);

fprintf('HPIPM returned with flag %d', return_flag);

if return_flag==0
    fprintf('-> QP solved! Solution:\n')
    qp_sol.print_C_struct()
else
    fprintf('-> Solver failed!')
end


% extract and print sol
fprintf('u =\n');
u = qp_sol.get_u();
for i=1:N+1
	u{i}
end

fprintf('x =\n');
for i=0:N
	x_tmp = qp_sol.get_x(i);
	x_tmp
end


% print solver statistics
fprintf('\nsolver statistics:\n\n');
fprintf('ipm return = %d\n\n', return_flag);
res_stat = solver.get_res_stat();
fprintf('ipm max res stat = %e\n\n', res_stat);
res_eq = solver.get_res_eq();
fprintf('ipm max res eq   = %e\n\n', res_eq);
res_ineq = solver.get_res_ineq();
fprintf('ipm max res ineq = %e\n\n', res_ineq);
res_comp = solver.get_res_comp();
fprintf('ipm max res comp = %e\n\n', res_comp);
iters = solver.get_iter();
fprintf('ipm iter = %d\n\n', iters);
stat = solver.get_stat();
fprintf('stat =\n');
fprintf('\talpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\n');
for ii=1:iters
	fprintf('\t%e\t%e\t%e\t%e\t%e\n', stat(ii,1), stat(ii,2), stat(ii,3), stat(ii,4), stat(ii,5));
end
fprintf('\n');

return

