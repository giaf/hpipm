% Matlab file to replace env.sh to allow compile_mex_all.m to run on Windows
% Len van Moorsel
% 2021.01.21

% set environment variables, that are originally set in env.sh

% set environment run variable to true to indicate this script was run
setenv('ENV_RUN','true')

% set paths to different folders used by HPIPM
setenv('HPIPM_MAIN_FOLDER',fileparts(fileparts(pwd))); % "$(pwd)/../.."
setenv('BLASFEO_MAIN_FOLDER',[fileparts(fileparts(fileparts(pwd))),'/blasfeo']); % "$(pwd)/../../../blasfeo"
setenv('MATLABPATH',[getenv('HPIPM_MAIN_FOLDER'),'/interfaces/matlab_octave/']); % $HPIPM_MAIN_FOLDER/interfaces/matlab_octave/
setenv('OCTAVE_PATH',[getenv('HPIPM_MAIN_FOLDER'),'/interfaces/matlab_octave/']); %$HPIPM_MAIN_FOLDER/interfaces/matlab_octave/
setenv('LD_LIBRARY_PATH',[getenv('HPIPM_MAIN_FOLDER'),'/lib/']); %$BLASFEO_MAIN_FOLDER/lib