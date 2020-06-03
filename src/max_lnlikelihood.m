function [optimum_logtheta, lnlikelihood_max] = max_lnlikelihood(samples,theta_searchbound)
% Function to find the maximum of the ln liklihood function

% This is done by minimising the minus of the ln likelihood function using
% a genetic algorithm search. This will give the optimum parameter theta
% for the Gaussian correlation function. Note theta is an array with length
% equal to the number of dimensions

%% Setup the GA search
% Create arrays of the theta search bounds sized as problem dimension
log_theta_lower = ones(1,samples.ndim).*theta_searchbound(1);
log_theta_upper = ones(1,samples.ndim).*theta_searchbound(2);
fprintf('\n----- Estimating Gaussian Hyperparameters -----\n');

% Use MATLAB's ga function to find the minimum of the negative ln likelihood
% function (therefore maximum ln likelihood).
% Need to create a nested function with only one input and one output. 

gaoptions = optimoptions('ga','UseParallel', true, 'UseVectorized', false, 'Display','iter',...
    'FunctionTolerance',1e-3, 'PopulationSize', 500, 'MaxGenerations', 2000);
[optimum_logtheta,min_neglnlikelihood] = ...
    ga(@lnlikelihood,samples.ndim,[],[],[],[],log_theta_lower,log_theta_upper,[],gaoptions);

% transpose to make it easier to see
optimum_logtheta = optimum_logtheta';
% max ln likelihood is the negative of the found minimum
lnlikelihood_max = -min_neglnlikelihood;

%% Nested function which calculates negative of ln likelihood for given theta
% GA is carried out on this function to find the theta which minimises the
% negative ln likelihood function. Nested function is used since MATLAB's
% ga only takes a function with one input and one output. In this way, the
% variable "sample" can be accessed and used in the function below.

    function neglnlike = lnlikelihood(log_theta)
        % Input= log10 of theta array
        % Output = neglnlike = Negative of the ln likelihood function which
        % is to be minimised
        
        o = ones(samples.npoint,1);
        z = zeros(samples.npoint * samples.ndim,1);
        one = [o;z];
        
        % Get correlation matrix between sample points
        % corrmat doesn't take log theta input so convert
        R_gek = corrmat(samples, 10.^log_theta);
        
        [U,q] = chol(R_gek); % Cholesky Factorisation
        
        if q~=0
            % R is ill conditioned so discard neglnlike in genetic algorighm search
            % by assigning a high value. This ensures we get a well conditioned R
            neglnlike = eps^-1;
        else            
            mu = (one'*(U\(U'\samples.output_aug)))/(one'*(U\(U'\one)));
            sighat=((samples.output_aug-one*mu)'*(U\(U'\(samples.output_aug-one*mu))))/samples.npoint_gek;
            
            % Sum lns of diagonal to find ln(det(R))
            LnDetR=2*sum(log(abs(diag(U))));
            % minus ln likelihood to minimise is:
            neglnlike = -0.5*(-samples.npoint_gek*log(sighat)-LnDetR); 
        end
        
    end

end