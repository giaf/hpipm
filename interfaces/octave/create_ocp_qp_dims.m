function dims = create_ocp_qp_dims( N )
% dims

% control horizon
dims.N = N;
% number of states
dims.nx = zeros(N+1, 1);
% number of inputs
dims.nu = zeros(N+1, 1);
% number of state box constraints
dims.nbx = zeros(N+1, 1);
% number of input box constraints
dims.nbu = zeros(N+1, 1);
% number of general constraints
dims.ng = zeros(N+1, 1);
% number of state box constraints which are softed
dims.nsbx = zeros(N+1, 1);
% number of input box constraints which are softed
dims.nsbu = zeros(N+1, 1);
% number of general constraints which are softed
dims.nsg = zeros(N+1, 1);

return;

