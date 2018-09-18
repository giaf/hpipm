% pyversion /usr/bin/python3
% py.sys.path

clear all
close all
clc

addpath('../../interfaces/matlab/hpipm_matlab')
import hpipm_matlab.*



% dims
N = 5;

tic
dims = hpipm_ocp_qp_dim(N);
tmp_time = toc
fprintf('create dim time %e\n', tmp_time);

tic
dims.set_nx([2, 2, 2, 2, 2, 2]);
tmp_time = toc
fprintf('set nx time %e\n', tmp_time);
dims.set_nu([1, 1, 1, 1, 1]);
dims.set_nbx(2, 0);
dims.set_nbx(2, 5);

dims.print_C_struct();



% data
%A = [1, 1; 0, 1];
A = [1, 0; 1, 1];
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
tmp_time = toc
fprintf('create qp time %e\n', tmp_time);

tic
qp.set_A({A, A, A, A, A});
tmp_time = toc
fprintf('create set A time %e\n', tmp_time);
qp.set_B({B, B, B, B, B});
%qp.set_b({b, b, b, b, b});

qp.set_Q({Q, Q, Q, Q, Q, Q});
qp.set_S({S, S, S, S, S});
qp.set_R({R, R, R, R, R});
qp.set_q({q, q, q, q, q, q});
qp.set_q(q, 0);
%qp.set_r({r, r, r, r, r});
qp.set_Jx(Jx, 0);
qp.set_lx(x0, 0);
qp.set_ux(x0, 0);
qp.set_Jx(Jx, 5);

qp.print_C_struct()



% qp sol
tic
qp_sol = hpipm_ocp_qp_sol(dims);
tmp_time = toc
fprintf('create qp_sol time %e\n', tmp_time);



% set up solver
tic
solver = hpipm_ocp_qp_solver(dims);
tmp_time = toc
fprintf('create solver time %e\n', tmp_time);



% solve qp
tic
%return_flag = solver.solve(qp, qp_sol);
solver.solve(qp, qp_sol);
tmp_time = toc
fprintf('solve time %e\n', tmp_time);


qp_sol.print_C_struct()



return



qp_data = hpipm_data();

A = [1, 0; 1, 1];
B = [0, 1].';
b = [0, 0].';

Q = [1, 0; 0, 1];
S = [0, 0].';
R = 1;
q = [1, 1].';
r = 0;

qp_data.A = {A, A, A, A, A};
qp_data.B = {B, B, B, B, B};
qp_data.b = {b, b, b, b, b};
qp_data.Q = {Q, Q, Q, Q, Q, Q};
qp_data.S = {S, S, S, S, S, S};
qp_data.R = {R, R, R, R, R, R};
qp_data.q = {q, q, q, q, q, q};
qp_data.r = {r, r, r, r, r, r};

x0 = [1, 1].';

qp_data.d_lb = {x0};
qp_data.d_ub = {x0};

qp_data.idxb = {[1, 2].'};

qp_dims = hpipm_dims();

qp_dims.nx   = [2, 2, 2, 2, 2, 2].';
qp_dims.nu   = [1, 1, 1, 1, 1, 0].';
qp_dims.nb   = [2, 0, 0, 0, 0, 0].';
qp_dims.nbx  = [2, 0, 0, 0, 0, 0].';
qp_dims.nbu  = [0, 0, 0, 0, 0, 0].';
qp_dims.ng   = [0, 0, 0, 0, 0, 0].';
qp_dims.ns   = [0, 0, 0, 0, 0, 0].';
qp_dims.nsbx = [0, 0, 0, 0, 0, 0].';
qp_dims.nsbu = [0, 0, 0, 0, 0, 0].';
qp_dims.nsg  = [0, 0, 0, 0, 0, 0].';
qp_dims.N    = 5;

% set up solver
solver = hpipm_solver(qp_dims, qp_data);

% solve qp
return_flag = solver.solve();

% print solution
solver.print_sol();

