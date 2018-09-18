classdef hpipm_ocp_qp_dim
	


	properties
		hp
		np
		py_dims
	end



	methods
		

		function obj = hpipm_ocp_qp_dim(N)
%			import hpipm_matlab.*;
			obj.np = py.importlib.import_module('numpy');
			obj.hp = py.importlib.import_module('hpipm_python');
			obj.py_dims = obj.hp.hpipm_ocp_qp_dim(int32(N));
		end


		function set_nx(obj, varargin)
			if(length(varargin)==1)
				nx = varargin{1};
				m = length(nx);
				py_nx = obj.np.zeros({int32(m), int32(1)}, obj.np.int32);
				for i=1:m
					py_nx.itemset(int32(i-1), int32(0), int32(nx(i)));
				end
				obj.py_dims.set_nx(py_nx);
			else
				nx = varargin{1};
				idx = varargin{2};
				obj.py_dims.set_nx(int32(nx), int32(idx));
			end
		end
	

		function set_nu(obj, varargin)
			if(length(varargin)==1)
				nu = varargin{1};
				m = length(nu);
				py_nu = obj.np.zeros({int32(m), int32(1)}, obj.np.int32);
				for i=1:m
					py_nu.itemset(int32(i-1), int32(0), int32(nu(i)));
				end
				obj.py_dims.set_nu(py_nu);
			else
				nu = varargin{1};
				idx = varargin{2};
				obj.py_dims.set_nu(int32(nu), int32(idx));
			end
		end
	

		function set_nbx(obj, varargin)
			if(length(varargin)==1)
				nbx = varargin{1};
				m = length(nbx);
				py_nbx = obj.np.zeros({int32(m), int32(1)}, obj.np.int32);
				for i=1:m
					py_nbx.itemset(int32(i-1), int32(0), int32(nbx(i)));
				end
				obj.py_dims.set_nbx(py_nbx);
			else
				nbx = varargin{1};
				idx = varargin{2};
				obj.py_dims.set_nbx(int32(nbx), int32(idx));
			end
		end
	

		function set_nbu(obj, varargin)
			if(length(varargin)==1)
				nbu = varargin{1};
				m = length(nbu);
				py_nbu = obj.np.zeros({int32(m), int32(1)}, obj.np.int32);
				for i=1:m
					py_nbu.itemset(int32(i-1), int32(0), int32(nbu(i)));
				end
				obj.py_dims.set_nbu(py_nbu);
			else
				nbu = varargin{1};
				idx = varargin{2};
				obj.py_dims.set_nbu(int32(nbu), int32(idx));
			end
		end
	

		function set_ng(obj, varargin)
			if(length(varargin)==1)
				ng = varargin{1};
				m = length(ng);
				py_ng = obj.np.zeros({int32(m), int32(1)}, obj.np.int32);
				for i=1:m
					py_ng.itemset(int32(i-1), int32(0), int32(ng(i)));
				end
				obj.py_dims.set_ng(py_ng);
			else
				ng = varargin{1};
				idx = varargin{2};
				obj.py_dims.set_ng(int32(ng), int32(idx));
			end
		end
	

		function set_ns(obj, varargin)
			if(length(varargin)==1)
				ns = varargin{1};
				m = length(ns);
				py_ns = obj.np.zeros({int32(m), int32(1)}, obj.np.int32);
				for i=1:m
					py_ns.itemset(int32(i-1), int32(0), int32(ns(i)));
				end
				obj.py_dims.set_ns(py_ns);
			else
				ns = varargin{1};
				idx = varargin{2};
				obj.py_dims.set_ns(int32(ns), int32(idx));
			end
		end
	

		function print_C_struct(obj)
			obj.py_dims.print_C_struct();
		end


	end

end
