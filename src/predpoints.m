function [pred] = predpoints(sample, param, options)
% Obtain prediction points for GEK
% input options for prediction and sample
% output prediction struct

% Create prediction points
% Need a unifrom distribution of points across the boundaries of each
% design parameter. Done using Halton sequences.

% Number of prediction points as specified by user
pred.npoint = options.npred;

% Halton sequence. Skip and Leap values defined here
skip = floor(rand*1e7);
% leap = nextprime(sample.ndim)*6; % integer multiple of next prime
leap = nthprime(sample.ndim+1)-1;
halton = haltonset(sample.ndim,'Skip',skip,'Leap',leap);
halton = scramble(halton,'RR2');

% Create prediction points
pred.raw = net(halton, pred.npoint);

% Map prediction points [0 1] to bounds of each parameter.
[pred.mapped] = map_samples(param, pred.raw, options);

if strcmp(options.objective,'batch')   
    % for batch options.objective, add original sample points to the prediction points
    % matrix. the GEK MSE at these points should be ~0    
    pred.mapped = vertcat(pred.mapped, sample.input);
    pred.npoint = pred.npoint + sample.npoint;

elseif strcmp(options.objective,'verify')
    % for verify options.objective, replace all SA values in prediction matrix with
    % nominal SA. Verify if GEK is producing same output for velocity options.objective
    % function than SU2 RANS at nomincal SA.
   
    % nominal SA values
    SAnom = [0.136,0.66667,0.622,0.41,0.3,2,7.1];
    SAnom = repmat(SAnom,pred.npoint,1);
    % populate prediction points
    pred.mapped = horzcat(SAnom,pred.mapped(:,param.x),pred.mapped(:,param.y));
end

end

