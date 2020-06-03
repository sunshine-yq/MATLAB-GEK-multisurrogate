function [mu, sighat] = kriging_mean(samples, R)
% Find kriging mean and sigma

%cholesky factorisation
[U] = chol(R);

% Kriging mean
o = ones(samples.npoint,1);
z = zeros(samples.npoint * samples.ndim,1);
one = [o;z];
mu = (one'*(U\(U'\samples.output_aug)))/(one'*(U\(U'\one)));
sighat=((samples.output_aug-one*mu)'*(U\(U'\(samples.output_aug-one*mu))))/samples.npoint_gek;

end

