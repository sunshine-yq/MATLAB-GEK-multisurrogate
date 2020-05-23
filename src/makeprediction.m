function [pred] = makeprediction(sample, pred, GEK)
% Make GEK predictions on points to find output and MSE

% Initialise arrays to store
predoutput = zeros(pred.npoint, 1);
predmse = zeros(pred.npoint, 1);

% de-struct for parfor loop
predpoint = pred.mapped;
theta = GEK.theta;
mu = GEK.mu;
sighat = GEK.sighat;

% Cholesky Factorisation
U = chol(GEK.R);

% one zero matrices
o = ones(sample.npoint,1);
z = zeros(sample.npoint * sample.ndim,1);
one = [o;z];

parfor ii = 1:pred.npoint
    
    % Find correlation between prediction point and samples
    r = corrmat_pred(sample, theta, predpoint(ii,:));
    
    % Make prediction and find error
    y_p = mu + r'*(U\(U'\(sample.output_aug-one*mu)));
    mse = sighat*(1-r'*(U\(U'\r))+((1-one'*(U\(U'\r)))^2/(one'*(U\(U'\one)))));
    
    % De-normalise and store
    predoutput(ii) = y_p*sample.output_sd + sample.output_mean;
    predmse(ii) = mse/sighat;
    
end

% MSE values smaller than eps assigned 0
predmse(predmse<eps)=0;

% put results back into main pred struct
pred.output = predoutput;
pred.mse = predmse;

% Sort prediction mse from max to min
[pred.mse_sortval, pred.mse_sortindex] = sort(pred.mse,'descend');

end

