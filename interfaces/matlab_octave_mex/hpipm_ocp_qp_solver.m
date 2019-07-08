classdef hpipm_ocp_qp_solver < handle
	
	properties
		C_dim
		C_arg
		C_ws
		time_ext
	end

	methods

		function obj = hpipm_ocp_qp_solver(dim, arg)
			obj.C_dim = dim.C_dim;
			obj.C_arg = arg.C_arg;
			% create struct in C
			obj.C_ws = ocp_qp_solver_create(obj.C_dim, obj.C_arg);
		end

		function solve(obj, qp, sol)
			obj.time_ext = ocp_qp_solver_solve(qp.C_qp, sol.C_sol, obj.C_arg, obj.C_ws);
		end

		function set(obj, field, value)
			if(strcmp(field, 'iter_max') || strcmp(field, 'tol_stat') || strcmp(field, 'tol_eq') || strcmp(field, 'tol_ineq') || strcmp(field, 'tol_comp') || strcmp(field, 'mu0') || strcmp(field, 'reg_prim'))
				ocp_qp_solver_arg_set(obj.C_arg, field, value);
			else
				disp(['hpipm_ocp_qp_solver.set: wrong field: ', field]);
				keyboard;
			end
		end

		function value = get(obj, field)
			if(strcmp(field, 'time_ext'))
				value = obj.time_ext;
			else
				value = ocp_qp_solver_get(obj.C_ws, field);
			end
		end

%		function print_C_struct(obj)
%			ocp_qp_print(obj.C_dim, obj.C_qp);
%		end

		function delete(obj)
			%disp('in destructor');
			ocp_qp_solver_destroy(obj.C_ws);
		end

	end
end



