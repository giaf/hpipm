function qp = create_ocp_qp( dims )
% qp

N = dims.N;
nx = dims.nx;
nu = dims.nu;
nbx = dims.nbx;
nbu = dims.nbu;
ng = dims.ng;
nsbx = dims.nsbx;
nsbu = dims.nsbu;
nsg = dims.nsg;

%
qp.A = {};
for ii=1:N
	qp.A{ii} = zeros(nx(ii+1), nx(ii));
end
%
qp.B = {};
for ii=1:N
	qp.B{ii} = zeros(nx(ii+1), nu(ii));
end
%
qp.b = {};
for ii=1:N
	qp.b{ii} = zeros(nx(ii+1), 1);
end
%
qp.Q = {};
for ii=1:N+1
	qp.Q{ii} = zeros(nx(ii), nx(ii));
end
%
qp.R = {};
for ii=1:N+1
	qp.R{ii} = zeros(nu(ii), nu(ii));
end
%
qp.S = {};
for ii=1:N+1
	qp.S{ii} = zeros(nu(ii), nx(ii));
end
%
qp.q = {};
for ii=1:N+1
	qp.q{ii} = zeros(nx(ii), 1);
end
%
qp.r = {};
for ii=1:N+1
	qp.r{ii} = zeros(nu(ii), 1);
end
%
qp.Jx = {};
for ii=1:N+1
	qp.Jx{ii} = zeros(nbx(ii), nx(ii));
end
%
qp.lx = {};
for ii=1:N+1
	qp.lx{ii} = zeros(nbx(ii), 1);
end
%
qp.ux = {};
for ii=1:N+1
	qp.ux{ii} = zeros(nbx(ii), 1);
end
%
qp.Ju = {};
for ii=1:N+1
	qp.Ju{ii} = zeros(nbu(ii), nu(ii));
end
%
qp.lu = {};
for ii=1:N+1
	qp.lu{ii} = zeros(nbu(ii), 1);
end
%
qp.uu = {};
for ii=1:N+1
	qp.uu{ii} = zeros(nbu(ii), 1);
end

return;
