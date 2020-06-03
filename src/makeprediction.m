function [predpoints] = makeprediction(samples, predpoints, GEK)
% Make GEK predictions on points to find output and MSE

% Initialise arrays to store
predoutput = zeros(predpoints.npoint, 1);
predmse = zeros(predpoints.npoint, 1);

% de-struct for parfor loop
pred = predpoints.mapped;
theta = GEK.theta;
mu = GEK.mu;
sighat = GEK.sighat;

% Cholesky Factorisation
U = chol(GEK.R);

% one zero matrices
o = ones(samples.npoint,1);
z = zeros(samples.npoint * samples.ndim,1);
one = [o;z];

parfor ii = 1:predpoints.npoint
    
    % Find correlation between prediction point and samples
    r = corrmat_pred(samples, theta, pred(ii,:));
    
    % Make prediction and find error
    y_p = mu + r'*(U\(U'\(samples.output_aug-one*mu)));
    mse = sighat*(1-r'*(U\(U'\r))+((1-one'*(U\(U'\r)))^2/(one'*(U\(U'\one)))));
    
    % De-normalise and store
    predoutput(ii) = y_p*samples.output_sd + samples.output_mean;
    predmse(ii) = mse/sighat;
    
end

% MSE values smaller than eps assigned 0
predmse(predmse<eps)=0;

% put results back into main struct
predpoints.output = predoutput;
predpoints.mse = predmse;

% Sort prediction mse from max to min
[predpoints.mse_sortval, predpoints.mse_sortindex] = sort(predpoints.mse,'descend');

end

