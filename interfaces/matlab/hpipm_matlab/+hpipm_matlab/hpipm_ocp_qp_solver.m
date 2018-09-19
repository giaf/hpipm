classdef hpipm_ocp_qp_solver
	


	properties
		hp
		np
		py_arg
		py_solver
	end



	methods
		

		function obj = hpipm_ocp_qp_solver(dim, arg)
%			import hpipm_matlab.*;
			obj.np = py.importlib.import_module('numpy');
			obj.hp = py.importlib.import_module('hpipm_python');
			obj.py_solver = obj.hp.hpipm_ocp_qp_solver(dim.py_dims, arg.py_arg);
			obj.py_arg = arg.py_arg;
		end


		function return_flag = solve(obj, qp, qp_sol)
			py_flag = obj.py_solver.solve(qp.py_qp, qp_sol.py_qp_sol);
			return_flag = int64(py_flag.real);
		end
	

	end

end
	
