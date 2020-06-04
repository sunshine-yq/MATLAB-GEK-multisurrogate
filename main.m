%% Gradient Enhanced Kriging Modelling
% Model the product of flow velocity magnitude and angle at XY for given
% coefficients of the SA turbulence model
% Generate model from SU2 simulations

% Amir Bagheri
% May 2020

clear; close all; addpath(genpath('./'));

%% Set Options for running the code

options.platform    = 'local'; % platform to run on (iridis/local)
options.iterno      = 2; %

options.writetofile  = true;
options.nnextsamples = 60;
options.nfiles       = 1;
options.theta        = 'theta01';
options.npredpoints  = 1000;
options.objective    = 'iterate';
options.xcludradmax  = 0.02; % maximum exclusion radius 
options.xcludtanh    = 3; % tanh factor p. larger = more space b/w samples

options.nsurrogate  = 10;
options.activesrrgt = 1;

%%check_options()

%% Run Program

% Prepare environment for run
prepare_env();

% Create parameter struct and set integers to each parameter
% This aslo sets the dimension number of the model
param.cb1 = 1; param.sig = 2; param.cb2 = 3; param.kar = 4;
param.cw2 = 5; param.cw3 = 6; param.cv1 = 7;

% Create baslines samples if at Iteration 1 and stop progressing
if options.iterno == 1
   baseline_samples(param, options);
   return;
end

% Read the samples
[samples] = read_io(param, options);

% Initialise the parallel run
init_parallel(options);

% Calculate the hyperparameters of the Gaussian Correlation Function
[GEK.theta, GEK.ln_likelihood] = hyper_param(samples, options);

% Find Correlation matrix and Kriging mean using the hyperparameters
[GEK.R] = corrmat(samples, GEK.theta);
[GEK.mu, GEK.sighat] = kriging_mean(samples, GEK.R);

% Generate prediction points for GEK prediction
[predictions] = generate_predpoints(samples, param, options);

% Make GEK predictions
fprintf('\n----- Making Predictions -----\n');
[predictions] = makeprediction(samples, predictions, GEK);
fprintf('-Complete\n');

% Find next iteration of sample points
if strcmp(options.objective, 'iterate')
    fprintf('\n+++++ Selecting Next Iteration Samples +++++\n');
    [nextsamples] = next_iteration(predictions, options);
    fprintf('-Complete\n');
end




