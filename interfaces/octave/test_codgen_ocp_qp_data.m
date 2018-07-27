%

% horizon length
N = 5;

% dims
dims = create_ocp_qp_dims(N);
% number of states
for ii=1:N+1
	dims.nx(ii) = 2;
end
% number of inputs
for ii=1:N
	dims.nu(ii) = 1;
end
% number of state box constraints
dims.nbx(1) = dims.nx(1);

%dims

% data
% 
A = [1 1; 0 1];
%
B = [0; 1];
%
b = [0; 0];
%
Q = [1 0; 0 1];
%
S = [0 0];
%
R = [1];
%
q = [1; 1];
%
r = [0];
%
x0 = [1; 1];
%
Jx0 = [1 0; 0 1];

qp = create_ocp_qp(dims);

%
for ii=1:N
	qp.A{ii} = A;
end
%
for ii=1:N
	qp.B{ii} = B;
end
%
for ii=1:N
	qp.b{ii} = b;
end
%
for ii=1:N+1
	qp.Q{ii} = Q;
end
%
for ii=1:N
	qp.S{ii} = S;
end
%
for ii=1:N
	qp.R{ii} = R;
end
%
for ii=1:N+1
	qp.q{ii} = q;
end
%
for ii=1:N
	qp.r{ii} = r;
end
%
qp.Jx{1} = Jx0;
%
qp.lx{1} = x0;
%
qp.ux{1} = x0;

%qp

codegen_ocp_qp_data(dims, qp);

