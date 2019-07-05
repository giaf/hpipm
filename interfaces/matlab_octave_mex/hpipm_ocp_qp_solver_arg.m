classdef hpipm_ocp_qp_solver_arg < handle
	
	properties
		dim
		C_dim
		C_arg
	end

	methods

		function obj = hpipm_ocp_qp_solver_arg(dim)
			obj.dim = dim;
			obj.C_dim = dim.C_dim;
			% create struct in C
			obj.C_arg = ocp_qp_solver_arg_create(obj.C_dim);
		end

		function set(obj, field, value)
			ocp_qp_solver_arg_set(obj.C_arg, field, value);
		end

%		function print_C_struct(obj)
%			ocp_qp_print(obj.C_dim, obj.C_qp);
%		end

%		function codegen(obj, file_name, mode)
%			ocp_qp_codegen(obj.C_dim, obj.C_qp, file_name, mode);
%		end

		function delete(obj)
			%disp('in destructor');
			ocp_qp_solver_arg_destroy(obj.C_arg);
		end

	end
end


