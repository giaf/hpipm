classdef hpipm_ocp_qp
	


	properties
		hp
		np
		py_qp
	end



	methods
		

		function obj = hpipm_ocp_qp(dims)
%			import hpipm_matlab.*;
			obj.np = py.importlib.import_module('numpy');
			obj.hp = py.importlib.import_module('hpipm_python');
			obj.py_qp = obj.hp.hpipm_ocp_qp(dims.py_dims);
		end
	

		function set_A(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					tmp_mat = py_mat.append(tmp_mat);
				end
				obj.py_qp.set_A(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_A(tmp_mat, int32(idx));
			end
		end
	

		function set_B(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_B(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_B(tmp_mat, int32(idx));
			end
		end
	

		function set_b(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_b(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_b(tmp_mat, int32(idx));
			end
		end
	

		function set_Q(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_Q(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_Q(tmp_mat, int32(idx));
			end
		end
	

		function set_S(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_S(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_S(tmp_mat, int32(idx));
			end
		end
	

		function set_R(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_R(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_R(tmp_mat, int32(idx));
			end
		end
	

		function set_q(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_q(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_q(tmp_mat, int32(idx));
			end
		end
	

		function set_r(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_r(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_r(tmp_mat, int32(idx));
			end
		end
	

		function set_Jx(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_Jx(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_Jx(tmp_mat, int32(idx));
			end
		end
	

		function set_lx(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					tmp_mat
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_lx(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_lx(tmp_mat, int32(idx));
			end
		end
	

		function set_ux(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_ux(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_ux(tmp_mat, int32(idx));
			end
		end
	

		function set_Ju(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_Ju(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_Ju(tmp_mat, int32(idx));
			end
		end
	

		function set_lu(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_lu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_lu(tmp_mat, int32(idx));
			end
		end
	

		function set_uu(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_uu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_uu(tmp_mat, int32(idx));
			end
		end
	

		function set_C(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_C(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_C(tmp_mat, int32(idx));
			end
		end
	

		function set_D(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_D(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_D(tmp_mat, int32(idx));
			end
		end
	

		function set_lg(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_lg(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_lg(tmp_mat, int32(idx));
			end
		end
	

		function set_ug(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_ug(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_ug(tmp_mat, int32(idx));
			end
		end
	

		function set_Zl(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_Zl(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_Zl(tmp_mat, int32(idx));
			end
		end
	

		function set_Zu(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_Zu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_Zu(tmp_mat, int32(idx));
			end
		end
	

		function set_zl(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_zl(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_zl(tmp_mat, int32(idx));
			end
		end
	

		function set_zu(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_zu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_zu(tmp_mat, int32(idx));
			end
		end
	

		function set_lls(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_lls(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_lls(tmp_mat, int32(idx));
			end
		end
	

		function set_lus(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_lus(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_lus(tmp_mat, int32(idx));
			end
		end
	

		function set_Jsu(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_Jsu(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_Jsu(tmp_mat, int32(idx));
			end
		end
	

		function set_Jsx(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_Jsx(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%				tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_Jsx(tmp_mat, int32(idx));
			end
		end
	

		function set_Jsg(obj, varargin)
			if(length(varargin)==1)
				mat = varargin{1};
				m = length(mat);
				py_mat = py.list({});
				for i=1:m
					mat0 = mat{i};
					n = size(mat0, 1);
					k = size(mat0, 2);
					tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
					for j=1:n
						for l=1:k
							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%							tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
						end
					end
%					tmp_mat = tmp_mat.transpose();
					py_mat.append(tmp_mat);
				end
				obj.py_qp.set_Jsg(py_mat);
			else
				mat0 = varargin{1};
				idx = varargin{2};
				n = size(mat0, 1);
				k = size(mat0, 2);
				tmp_mat = obj.np.zeros({int32(n), int32(k)}, obj.np.double);
				for j=1:n
					for l=1:k
						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(j,l)));
%						tmp_mat.itemset(int32(j-1), int32(l-1), double(mat0(l,j)));
					end
				end
%					tmp_mat = tmp_mat.transpose();
				obj.py_qp.set_Jsg(tmp_mat, int32(idx));
			end
		end
	

		function print_C_struct(obj)
			obj.py_qp.print_C_struct();
		end


	end

end



