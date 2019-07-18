classdef hpipm_ocp_qp < handle
	
	properties
		dim
		C_dim
		C_qp
	end

	methods

		function obj = hpipm_ocp_qp(dim)
			obj.dim = dim;
			obj.C_dim = dim.C_dim;
			% create struct in C
			obj.C_qp = ocp_qp_create(obj.C_dim);
		end

		function set(varargin)
			if nargin==4
				obj = varargin{1};
				field = varargin{2};
				value = varargin{3};
				stage0 = varargin{4};
				ocp_qp_set(obj.C_qp, field, value, stage0);
			elseif nargin==5
				obj = varargin{1};
				field = varargin{2};
				value = varargin{3};
				stage0 = varargin{4};
				stage1 = varargin{5};
				ocp_qp_set(obj.C_qp, field, value, stage0, stage1);
			else
				disp('hpipm_ocp_qp.set: wrong number of input arguments (3 or 4 allowed)');
			end
		end

		function print_C_struct(obj)
			ocp_qp_print(obj.C_dim, obj.C_qp);
		end

		function codegen(obj, file_name, mode)
			ocp_qp_codegen(obj.C_dim, obj.C_qp, file_name, mode);
		end

		function delete(obj)
			%disp('in destructor');
			ocp_qp_destroy(obj.C_qp);
		end

	end
end

