function compile_mex_ocp_qp()

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	disp('ERROR: env.sh has not been sourced! Before executing this example, run:');
	disp('source env.sh');
	return;
end

% get acados folder
hpipm_folder = getenv('HPIPM_MAIN_FOLDER');
blasfeo_folder = getenv('BLASFEO_MAIN_FOLDER');
mex_flags = getenv('HPIPM_MEX_FLAGS');

% set paths
hpipm_mex_folder = [hpipm_folder, '/interfaces/matlab_octave/'];
hpipm_include = ['-I', hpipm_folder, '/include'];
hpipm_lib = ['-L', hpipm_folder, '/lib'];
blasfeo_include = ['-I', blasfeo_folder, '/include'];
blasfeo_lib = ['-L', blasfeo_folder, '/lib'];

% compile mex
mex_files = { ...
	[hpipm_mex_folder, 'ocp_qp_dim_create.c'], ...
	[hpipm_mex_folder, 'ocp_qp_dim_destroy.c'], ...
	[hpipm_mex_folder, 'ocp_qp_dim_set.c'], ...
	[hpipm_mex_folder, 'ocp_qp_dim_print.c'], ...
	[hpipm_mex_folder, 'ocp_qp_dim_codegen.c'], ...
	[hpipm_mex_folder, 'ocp_qp_create.c'], ...
	[hpipm_mex_folder, 'ocp_qp_destroy.c'], ...
	[hpipm_mex_folder, 'ocp_qp_set.c'], ...
	[hpipm_mex_folder, 'ocp_qp_print.c'], ...
	[hpipm_mex_folder, 'ocp_qp_codegen.c'], ...
	[hpipm_mex_folder, 'ocp_qp_sol_create.c'], ...
	[hpipm_mex_folder, 'ocp_qp_sol_destroy.c'], ...
	[hpipm_mex_folder, 'ocp_qp_sol_get.c'], ...
	[hpipm_mex_folder, 'ocp_qp_sol_print.c'], ...
	[hpipm_mex_folder, 'ocp_qp_solver_arg_create.c'], ...
	[hpipm_mex_folder, 'ocp_qp_solver_arg_destroy.c'], ...
	[hpipm_mex_folder, 'ocp_qp_solver_arg_set.c'], ...
	[hpipm_mex_folder, 'ocp_qp_solver_arg_codegen.c'], ...
	[hpipm_mex_folder, 'ocp_qp_solver_create.c'], ...
	[hpipm_mex_folder, 'ocp_qp_solver_destroy.c'], ...
	[hpipm_mex_folder, 'ocp_qp_solver_solve.c'], ...
	[hpipm_mex_folder, 'ocp_qp_solver_get.c'], ...
	} ;

if is_octave()
	if exist('cflags_octave.txt')==0
		diary 'cflags_octave.txt'
		diary on
		mkoctfile -p CFLAGS
		diary off
		input_file = fopen('cflags_octave.txt', 'r');
		cflags_tmp = fscanf(input_file, '%[^\n]s');
		fclose(input_file);
		cflags_tmp = [cflags_tmp, ' -std=c99'];
		input_file = fopen('cflags_octave.txt', 'w');
		fprintf(input_file, '%s', cflags_tmp);
		fclose(input_file);
	end
	input_file = fopen('cflags_octave.txt', 'r');
	cflags_tmp = fscanf(input_file, '%[^\n]s');
	fclose(input_file);
	setenv('CFLAGS', cflags_tmp);
end

for ii=1:length(mex_files)
	disp(['compiling ', mex_files{ii}])
	if is_octave()
%		mkoctfile -p CFLAGS
		mex(hpipm_include, blasfeo_include, hpipm_lib, blasfeo_lib, '-lhpipm', '-lblasfeo', mex_files{ii});
	else
		mex(mex_flags, 'CFLAGS=\$CFLAGS -std=c99', hpipm_include, blasfeo_include, hpipm_lib, blasfeo_lib, '-lhpipm', '-lblasfeo', mex_files{ii});
	end
end
