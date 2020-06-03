function [r_gek] = corrmat_pred(samples, theta, pred)
% Correlation matrix between prediction point and sample points

%% Initialise

r_gek = zeros(samples.npoint_gek,1);
p = 2; % Power of RBF. 2 = Gaussian Correlation

%% Standard kriging correlation matrix
r = ones(samples.npoint,1);
for z=1:samples.ndim
    r_int = zeros(samples.npoint,1);
    for i = 1:samples.npoint
        d = abs(samples.input(i,z) - pred(z));
        r_int(i) = exp(-theta(z)*d^p);
    end
    r = r.*r_int; % Multiply r_int successively for all dimensions
end

% Populate r_gek
for i=1:samples.npoint
    r_gek(i) = r(i);
end

%% Derivative correlation matrices
for z = 1:samples.ndim
    r = zeros(samples.npoint,1);
    for i = 1:samples.npoint
        d = samples.input(i,z) - pred(z);
        r(i) = -p*theta(z)*d*r_gek(i);
    end
    
    % Populate r_gek
    for i=1:samples.npoint
        r_gek(i+z*samples.npoint) = r(i);
    end
end

end

