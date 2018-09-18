classdef hpipm_ocp_qp_sol
	


	properties
		hp
		np
		py_qp_sol
	end



	methods
		

		function obj = hpipm_ocp_qp_sol(dims)
%			import hpipm_matlab.*;
			obj.np = py.importlib.import_module('numpy');
			obj.hp = py.importlib.import_module('hpipm_python');
			obj.py_qp_sol = obj.hp.hpipm_ocp_qp_sol(dims.py_dims);
		end
	

		function print_C_struct(obj)
			obj.py_qp_sol.print_C_struct();
		end


	end

end
