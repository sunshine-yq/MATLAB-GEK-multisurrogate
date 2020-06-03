function [R_gek] = corrmat(samples, theta)
% Find the correlation between the sample points using Gausian correlation function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R is made of multiple sub-matrices. Each sub-matrix is constructed
% seperately and then are all amalgamated in the end for speed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialise
n = samples.npoint;
n_gek = samples.npoint_gek;
R_gek = zeros(n_gek);

p = 2; % Power of RBF. 2 = Gaussian Correlation. Fixed

%% The standard kriging correlation matrix for all samples

R = ones(n,n);
for z=1:samples.ndim
    R_int = zeros(n,n); % Intermediary Correlation function for each dimension
    x = samples.input(:,z);
    k = 2; % Only populate upper triangle as matrix is symmetric
    for i=1:n-1
        for j=k:n
            d = abs(x(i) - x(j));
            R_int(i,j) = exp(-theta(z)*d^p);
        end
        k = k+1;
        % Another way same thing more efficient
%         d = abs(x(i) - x(k:end));
%         R_int(i,k:end) = exp(-theta(z).*d.^p)';
%         k = k+1;        
    end
    R_int = R_int + R_int' + eye(n); % Add transpose to symmetric matrix
    
    R = R.*R_int; % Multiply by R_int successively for all dimensions
end

% Add small diagonal to make matrix less ill-conditioned
R = R + eye(n)*eps;

% Populate R_gek
for i=1:n
    for j=1:n
        R_gek(i,j) = R(i,j);
    end
end

%% The first derivative correlation matrices

for z = 1:samples.ndim
    R = zeros(n);
    k = 2; % Only populate upper triangle as matrix is symmetric
    for i = 1:n-1
        for j = k:n
            R(i,j) = +p*theta(z)*(samples.input(i,z)-samples.input(j,z)) ...
                * R_gek(i,j);
        end
        k = k+1;
    end
    % Because the derivative correlation sub-matrix doensn't have abs(d) in
    % its formulation, it is symmetric but with opposite sign for the
    % transponsed arrays. Therefore we add negative transpose.
    R = R - R';
    
    % Populate R_gek with each matrix and its transpose (which is just the
    % same matrix with the opposite sign, since it was symmetric to begin
    % with)
    for i=1:n
        for j=1:n
            R_gek(i, j+z*n) = R(i,j);
            R_gek(i+z*n, j) = -R(i,j);
        end
    end
end

%% The second derivative matrices

for z1=1:samples.ndim
    for z2=1:samples.ndim
        R = zeros(n);
        k = 1; % Only populate upper triangle as matrix is symmetric
        for i=1:n
            for j=k:n
                if z1==z2
                    R(i,j) = ...
                        p*theta(z1)*(1-p*theta(z1)*((samples.input(i,z1)-samples.input(j,z2))^p))...
                        *R_gek(i,j);
                else
                    R(i,j) = ...
                        -p*p*theta(z1)*theta(z2)*(samples.input(i,z1)-samples.input(j,z1))...
                        *(samples.input(i,z2)-samples.input(j,z2))*R_gek(i,j);
                end
            end
            k = k+1;
        end
        % Add transpose to create symmetric matrix
        % Since diagonal elements are also calculated, subtract them to
        % avoid adding twice. diag is called twice to create diagonal
        % matrix of diagonal elements
        R = R + R' - diag(diag(R));
        
        % Populate R_gek with each matrix
        for i=1:n
            for j=1:n
                R_gek(i+z1*n, j+z2*n) = R(i,j);
            end
        end
    end
end
% Discard computationally small values for numerical stability
R_gek(abs(R_gek)<eps) = 0;
end



