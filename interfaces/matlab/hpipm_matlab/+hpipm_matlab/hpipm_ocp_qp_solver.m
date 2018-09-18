classdef hpipm_ocp_qp_solver
	


	properties
		hp
		np
		py_qp_solver
	end



	methods
		

		function obj = hpipm_ocp_qp_solver(dims)
%			import hpipm_matlab.*;
			obj.np = py.importlib.import_module('numpy');
			obj.hp = py.importlib.import_module('hpipm_python');
			obj.py_qp_solver = obj.hp.hpipm_ocp_qp_solver(dims.py_dims);
		end


		function solve(obj, qp, qp_sol)
			obj.py_qp_solver.solve(qp.py_qp, qp_sol.py_qp_sol);
		end
	

	end

end
	
