%% Gradient Enhanced Kriging Modelling
% Model the product of flow velocity magnitude and angle at XY for given
% coefficients of the SA turbulence model
% Generate model from SU2 simulations

% Amir Bagheri
% Feb 2020

clear; close all; rng('shuffle');

% Set path
addpath(genpath('./'));

%% Set Options for running the code

% General 
options.platform  = 'local'; % platform to run on (iridis/local)
options.nfiles    = 9; % Number of files to read from samples folder
options.theta     = 'theta11'; % theta file. If left blank found using GA
options.objective = 'batch'; % New sample "batch" or "verify" existing GEK prediction
options.npred     = 1000; % number of prediction points

% Global XY boundaries
options.globalx = [0.76 1.1];
options.globaly = [0 0.04146];

% Next sample batch
options.nbatch      = 200; % number of next sample batch points
options.batchmaxrad = 0.02; % maximum exclusion radius 
options.batchtanh   = 3; % tanh factor p. larger = more space b/w samples
options.batchxbound = []; % xy bounds to reduce window size
options.batchybound = [];
options.writebatch  = false; % Write next sample batch to file

%% Run program

% Check options and set defaults if unspecified by user
options = check_options(options);

% Initialise the parallel run
init_parallel(options);

% Read SU2 input and output from samples folder
[sample, param] = read_io(options);
if sample.nfiles == 0, return; end % terminate here if only generating baseline samples

% Calculate the hyperparameters theta of the Gaussian Correlation Function
tic;
[GEK.theta, GEK.ln_likelihood] = hyper_param(sample, options);
time.hyper = toc/60;

% Find Correlation matrix and Kriging mean using the hyperparameters
[GEK.R] = corrmat(sample, GEK.theta);
[GEK.mu, GEK.sighat] = kriging_mean(sample, GEK.R);

% Generate prediction points at which to predict GEK output
[pred] = predpoints(sample, param, options);

% Make GEK predictions
fprintf('\n----- Making Predictions -----\n');
tic;
[pred] = makeprediction(sample, pred, GEK);
time.prediction = toc/60;
fprintf('-Complete\n');

% Find next batch of sample points
if strcmp(options.objective, 'batch')
    fprintf('\n+++++ Selecting Batch +++++\n');
    [batch, pool] = nextbatch(sample, pred, param, options);
    fprintf('-Complete\n');
else
    batch = []; pool = [];
end

% Generate plots if on local
if strcmp(options.platform, 'local')
    plotgek(sample, param, pred, batch, pool, options)
end

% Save the workspace variables if on Iridis
if strcmp(options.platform, 'iridis')
    save(sprintf('../iridisout/allvars_%s',options.objective));
end

% Confirm success
fprintf('\n***** All Complete *****\n');
