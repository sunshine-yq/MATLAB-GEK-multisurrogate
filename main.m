%% Gradient Enhanced Kriging Modelling
% Model the product of flow velocity magnitude and angle at XY for given
% coefficients of the SA turbulence model
% Generate model from SU2 simulations

% Amir Bagheri
% May 2020

clear; close all; rng('shuffle');

% Set path
addpath(genpath('./'));

%% Set Options for running the code

options.platform  = 'local'; % platform to run on (iridis/local)
options.iterno = 0; %
options.nmodel = 10;
options.writetofile = true;
options.nsamplesnew = 60;

%%check_options()

%% Run Program

% Create parameter struct and set integers to each parameter
param.cb1 = 1; param.sig = 2; param.cb2 = 3; param.kar = 4;
param.cw2 = 5; param.cw3 = 6; param.cv1 = 7;

% Create baslines samples if at Iteration 0
if options.iterno == 0
   [samples] = baseline_samples(param, options); 
end



% Initialise the parallel run
% init_parallel(options);


