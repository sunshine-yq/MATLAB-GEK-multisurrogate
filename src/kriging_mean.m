function [mu, sighat] = kriging_mean(sample, R)
%find kriging mean and sigma
% input sample struct and correlation matrix

%cholesky factorisation
[U] = chol(R);

% Kriging mean
o = ones(sample.npoint,1);
z = zeros(sample.npoint * sample.ndim,1);
one = [o;z];
mu = (one'*(U\(U'\sample.output_aug)))/(one'*(U\(U'\one)));
sighat=((sample.output_aug-one*mu)'*(U\(U'\(sample.output_aug-one*mu))))/sample.npoint_gek;

end

