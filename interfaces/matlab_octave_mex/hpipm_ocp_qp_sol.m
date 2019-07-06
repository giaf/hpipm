classdef hpipm_ocp_qp_sol < handle
	
	properties
		dim
		C_dim
		C_sol
	end

	methods

		function obj = hpipm_ocp_qp_sol(dim)
			obj.dim = dim;
			obj.C_dim = dim.C_dim;
			% create struct in C
			obj.C_sol = ocp_qp_sol_create(obj.C_dim);
		end

		function value = get(varargin)
			if nargin==3
				obj = varargin{1};
				field = varargin{2};
				stage0 = varargin{3};
				value = ocp_qp_sol_get(obj.C_sol, field, stage0);
			elseif nargin==4
				obj = varargin{1};
				field = varargin{2};
				stage0 = varargin{3};
				stage1 = varargin{4};
				value = ocp_qp_sol_get(obj.C_sol, field, stage0, stage1);
			else
				disp('hpipm_ocp_qp_sol.get: wrong number of input arguments (2 or 3 allowed)');
				keyboard;
			end
		end

		function print_C_struct(obj)
			ocp_qp_sol_print(obj.C_dim, obj.C_sol);
		end

%		function codegen(obj, file_name, mode)
%			ocp_qp_codegen(obj.C_dim, obj.C_qp, file_name, mode);
%		end

		function delete(obj)
			%disp('in destructor');
			ocp_qp_sol_destroy(obj.C_sol);
		end

	end
end


