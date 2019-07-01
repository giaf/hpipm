classdef hpipm_ocp_qp_dim < handle
	
	properties
		C_dim
	end

	methods

		function obj = hpipm_ocp_qp_dim(N)
			% create dims struct in C
			obj.C_dim = ocp_qp_dim_create(N);
		end

		function set(varargin)
			if nargin==4
				obj = varargin{1};
				field = varargin{2};
				value = varargin{3};
				stage0 = varargin{4};
				ocp_qp_dim_set(obj.C_dim, field, value, stage0);
			elseif nargin==5
				obj = varargin{1};
				field = varargin{2};
				value = varargin{3};
				stage0 = varargin{4};
				stage1 = varargin{5};
				ocp_qp_dim_set(obj.C_dim, field, value, stage0, stage1);
			else
				disp('hpipm_ocp_qp_dim.set: wrong number of input arguments (3 or 4 allowed)');
			end
		end

		function print_C_struct(obj)
			ocp_qp_dim_print(obj.C_dim);
		end

		function codegen(obj, file_name, mode)
			ocp_qp_dim_codegen(obj.C_dim, file_name, mode);
		end

		function delete(obj)
			%disp('in destructor');
			ocp_qp_dim_destroy(obj.C_dim);
		end

	end
end
