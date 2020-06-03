function [] = prepare_env(~)
% Prepare conditions for run

rng('shuffle');
warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');

% Change the current folder to the folder of main.m.
if(~isdeployed)
  cd(fileparts(which('main.m')));
end
end

