%% Gradient Enhanced Kriging Modelling
% Model the product of flow velocity magnitude and angle at XY for given
% coefficients of the SA turbulence model
% Generate model from SU2 simulations

% Amir Bagheri
% May 2020

clear; close all; addpath(genpath('./'));

%% Set Options for running the code

options.nsurrogates  = 10;
options.activesrrgt  = 1;

options.platform     = 'local';
options.objective    = 'iterate';

options.nfiles       = 1;
options.npredpoints  = 1000;
options.nnextsamples = 100;
options.theta        = 'theta01';

options.writetofile  = true;

check_options(options);

%% Run Program

% Prepare environment for run
prepare_env();

% Create parameter struct and set integers to each parameter
% This aslo sets the dimension number of the model
param.cb1 = 1; param.sig = 2; param.cb2 = 3; param.kar = 4;
param.cw2 = 5; param.cw3 = 6; param.cv1 = 7;

% Create baslines samples if at first iteration and stop progressing
if options.nfiles == 0
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
[predictions] = make_prediction(samples, predictions, GEK);
fprintf('-Complete\n');

% Find next iteration of sample points
if strcmp(options.objective, 'iterate')
    fprintf('\n+++++ Selecting Next Iteration Samples +++++\n');
    [nextsamples] = next_iteration(predictions, options);
    fprintf('-Complete\n');
end

% Generate plots if on local
if strcmp(options.platform, 'local')
    plotgek(samples, param, predictions, nextsamples, options)
end

% Save the workspace variables if on Iridis
if strcmp(options.platform, 'iridis')
    save(sprintf('Iridisout/allvars_%s',options.objective));
end

% Confirm success
fprintf('\n***** All Complete *****\n');




