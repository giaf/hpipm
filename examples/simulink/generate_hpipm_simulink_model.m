% create and open the model
modelname = 'hpipm_simulink_getting_started_gen';
hpipm_solver_name = 'hpipm_solver_sfunction';

% create and open the model
open_system(new_system(modelname));

set_param(modelname,'PreloadFcn','x0 = [3.5 3.5 0 0 0 0 0 0]; nu = 3; nx = 8; Ts = 0.1;')

% add HPIPM solver block
add_block('simulink/User-Defined Functions/S-Function', [modelname,'/hpipm_solver']);
set_param([modelname,'/hpipm_solver'],'FunctionName', hpipm_solver_name);
set_param([modelname,'/hpipm_solver'],'position', [50    50    300    200]);


% add unit delay
add_block('simulink/Discrete/Unit Delay', [modelname,'/delay']);
set_param([modelname,'/delay'],'InitialCondition', 'x0');
set_param([modelname,'/delay'],'SampleTime', 'Ts');
set_param([modelname,'/delay'],'position', [150    300    200    350]);


% add selectors
add_block('simulink/Signal Routing/Selector', [modelname,'/S1']);
set_param([modelname,'/S1'],'position', [350    50    375    75]);
set_param([modelname,'/S1'],'InputPortWidth', '90');


add_block('simulink/Signal Routing/Selector', [modelname,'/S2']);
set_param([modelname,'/S2'],'position', [350    100    375    125]);
set_param([modelname,'/S2'],'InputPortWidth', '248');

% add scope
add_block('simulink/Commonly Used Blocks/Scope', [modelname,'/trajectories scope']);
set_param([modelname,'/trajectories scope'],'position', [500    50    525    75]);
set_param([modelname,'/trajectories scope'],'NumInputPorts', '2');


add_block('simulink/Commonly Used Blocks/Scope', [modelname,'/info scope']);
set_param([modelname,'/info scope'],'position', [500    100    525    125]);
set_param([modelname,'/info scope'],'NumInputPorts', '2');

% Simulink.BlockDiagram.arrangeSystem 
% Simulink.BlockDiagram.routeLine

% connect blocks
add_line(modelname,'S1/1','trajectories scope/1');
add_line(modelname,'S2/1','trajectories scope/2');
add_line(modelname,'hpipm_solver/1','S1/1');
add_line(modelname,'hpipm_solver/2','S2/1');
add_line(modelname,'hpipm_solver/3','info scope/1');
add_line(modelname,'hpipm_solver/4','info scope/2');
add_line(modelname,'S2/1','delay/1');
add_line(modelname,'delay/1','hpipm_solver/1');

% save the model
save_system(modelname);