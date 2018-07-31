% pyversion /usr/bin/python3
% py.sys.path

hp = py.importlib.import_module('hpipm_python');
np = py.importlib.import_module('numpy');

qp_data = hp.hpipm_data();

A = m2py([1, 0; 1, 1].', np);
B = m2py([0, 1].', np);
b = m2py([0, 0].', np);

Q = m2py([1, 0; 0, 1].', np);
S = m2py([0, 0].', np);
R = m2py([1].', np);
q = m2py([1, 1].', np);
r = m2py([0].', np);

qp_data.A = py.list({A, A, A, A, A});
qp_data.B = py.list({B, B, B, B, B});
qp_data.b = py.list({b, b, b, b, b});
qp_data.Q = py.list({Q, Q, Q, Q, Q, Q});
qp_data.S = py.list({S, S, S, S, S, S});
qp_data.R = py.list({R, R, R, R, R, R});
qp_data.q = py.list({q, q, q, q, q, q});
qp_data.r = py.list({r, r, r, r, r, r});

x0 = m2py([1, 1].', np);

qp_data.d_lb = py.list({x0});
qp_data.d_ub = py.list({x0});

qp_data.idxb = py.list({m2py([1, 2].', np)});

qp_dims = hp.hpipm_dims();

qp_dims.nx   = m2py([2, 2, 2, 2, 2, 2].', np);
qp_dims.nu   = m2py([1, 1, 1, 1, 1, 0].', np);
qp_dims.nb   = m2py([2, 0, 0, 0, 0, 0].', np);
qp_dims.nbx  = m2py([2, 0, 0, 0, 0, 0].', np);
qp_dims.nbu  = m2py([0, 0, 0, 0, 0, 0].', np);
qp_dims.ng   = m2py([0, 0, 0, 0, 0, 0].', np);
qp_dims.ns   = m2py([0, 0, 0, 0, 0, 0].', np);
qp_dims.nsbx = m2py([0, 0, 0, 0, 0, 0].', np);
qp_dims.nsbu = m2py([0, 0, 0, 0, 0, 0].', np);
qp_dims.nsg  = m2py([0, 0, 0, 0, 0, 0].', np);
qp_dims.N    = int32(5);

% set up solver
solver = hp.hpipm_solver(qp_dims, qp_data);

% solve qp
return_flag = solver.solve();

solver.print_sol();

