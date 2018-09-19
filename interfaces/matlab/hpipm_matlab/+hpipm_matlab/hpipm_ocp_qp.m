classdef hpipm_ocp_qp
	


	properties
		hp
		np
		py_qp
	end



	methods
		

		function obj = hpipm_ocp_qp(dims)
			obj.np = py.importlib.import_module('numpy');
			obj.hp = py.importlib.import_module('hpipm_python');
			obj.py_qp = obj.hp.hpipm_ocp_qp(dims.py_dims);
		end
	

		function set_A(obj, varargin)
			import hpipm_matlab.*
			if (length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_A(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_A(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_B(obj, varargin)
			import hpipm_matlab.*
			if (length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_B(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_B(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_b(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_b(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_b(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_Q(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_Q(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_Q(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_S(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_S(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_S(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_R(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_R(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_R(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_q(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_q(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_q(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_r(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_r(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_r(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_Jx(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_Jx(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_Jx(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_lx(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_lx(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_lx(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_ux(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_ux(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_ux(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_Ju(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_Ju(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_Ju(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_lu(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_lu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_lu(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_uu(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_uu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_uu(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_C(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_C(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_C(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_D(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_D(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_D(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_lg(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_lg(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_lg(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_ug(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_ug(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_ug(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_Zl(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_Zl(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_Zl(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_Zu(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_Zu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_Zu(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_zl(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_zl(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_zl(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_zu(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_zu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_zu(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_lls(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_lls(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_lls(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_lus(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_lus(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_lus(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_Jsu(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_Jsu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_Jsu(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_Jsx(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_Jsx(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_Jsx(m2py(mat0, obj.np), int32(idx));
			end
		end


		function set_Jsg(obj, varargin)
			import hpipm_matlab.*
			if(length(varargin)==1)
				mat = varargin{1};
				py_mat = py.list({});
				for i=1:length(mat)
					mat0 = mat{i};
					py_mat.append(m2py(mat0, obj.np));
				end
				obj.py_qp.set_Jsg(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				obj.py_qp.set_Jsg(m2py(mat0, obj.np), int32(idx));
			end
		end


		function print_C_struct(obj)
			obj.py_qp.print_C_struct();
		end


	end

end



