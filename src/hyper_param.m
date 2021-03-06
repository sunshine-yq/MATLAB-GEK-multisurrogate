function [theta, ln_likelihood] = hyper_param(samples, options)
% Obtain Gaussian Hyperparameter theta
% We assume Gaussian, therefore only theta is unknown and p=2.

% Check option whether to load theta or Genetic Algorithm
if ~isempty(options.theta)
    
    fprintf('\nLoading theta file %s \n', options.theta);
    theta=load(fullfile(sprintf('Optimum_Theta/M%.2i',options.activesrrgt),options.theta));    
    theta = theta.theta;
    % find ln linkelihood for this value of theta
    R = corrmat(samples, theta); % Get correlation matrix
    [ln_likelihood] = get_lnlikelihood(R, samples);
    
else
    % Do Genetic Algorithm search for best theta to maximise ln likelihood
    % Define lower and upper search boundary of the log10(theta)
    theta_searchbound = [-6,6];
    
    % Find optimum log10(theta) which maximises the ln likelihood function
    % theta is an array with same length as number of dimension; a theta is
    % found for each dimension.
    [optimum_logtheta, ln_likelihood] = max_lnlikelihood(samples,theta_searchbound);

    % Convert log10 to real
    theta = 10.^(optimum_logtheta);   
    
end

% Print theta to screen
% list of parameters in order
paras = ["cb1";"sig";"cb2";"kar";"cw2";"cw3";"cv1"];
[opthe, oppar]=sort(theta,'descend');
fprintf('\n');
disp('--------- Theta ---------')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Parameter    Theta(i)')
for i=1:samples.ndim
    line = sprintf('%s          %.3e\n',paras(oppar(i)), opthe(i));
    fprintf(line)
end
disp('~~~~~~~~~~~~~~~~~~~~~~~~~')
fprintf(sprintf('ln Likelihood = %.3e\n',ln_likelihood))
end

