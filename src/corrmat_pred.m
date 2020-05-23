function [r_gek] = corrmat_pred(sample, theta, pred)
% Correlation matrix between prediction point and sample points

% == Inputs ==
% sample = struct of sample points
% theta = Gaussian correlation function parameter vector. one theta for
% each dimension.
% pred = point to be predicted. array of size ndim (point in each dim)

% == Outputs ==
% r_gek = N(d+1)x1 Correlation matrix between samlple points and
% prediction point.

%% Initialise

r_gek = zeros(sample.npoint_gek,1);
p = 2; % Power of RBF. 2 = Gaussian Correlation

%% Standard kriging correlation matrix
r = ones(sample.npoint,1);
for z=1:sample.ndim
    r_int = zeros(sample.npoint,1);
    for i = 1:sample.npoint
        d = abs(sample.input(i,z) - pred(z));
        r_int(i) = exp(-theta(z)*d^p);
    end
    r = r.*r_int; % Multiply r_int successively for all dimensions
end

% Populate r_gek
for i=1:sample.npoint
    r_gek(i) = r(i);
end

%% Derivative correlation matrices
for z = 1:sample.ndim
    r = zeros(sample.npoint,1);
    for i = 1:sample.npoint
        d = sample.input(i,z) - pred(z);
        r(i) = -p*theta(z)*d*r_gek(i);
    end
    
    % Populate r_gek
    for i=1:sample.npoint
        r_gek(i+z*sample.npoint) = r(i);
    end
end

end

