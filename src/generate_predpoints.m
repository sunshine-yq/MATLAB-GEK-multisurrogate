function [predictions] = generate_predpoints(samples, param, options)
% Generate prediction points for GEK

% Need a unifrom distribution of points across the boundaries of each
% design parameter. Done using Halton sequences.

% Number of prediction points as specified by user
predictions.npoint = options.npredpoints;

% Halton sequence. Skip and Leap values defined here
skip = floor(rand*1e7);
leap = nthprime(samples.ndim+1)-1;
halton = haltonset(samples.ndim,'Skip',skip,'Leap',leap);
halton = scramble(halton,'RR2');

% Create prediction points between [0,1] (used in GEK prediction)
predictions.raw = net(halton, predictions.npoint);

% Add original sample points to the prediction points
% matrix. GEK MSE at these points should be ~0
predictions.raw = vertcat(predictions.raw, samples.input);
predictions.npoint = predictions.npoint + samples.npoint;

% Map prediction points to bounds of each parameter (used only in plot)
[predictions.mapped] = map_samples(param, predictions.raw);

end

