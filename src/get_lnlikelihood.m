function [lnlikelihood] = get_lnlikelihood(R, samples)
% Function to get the lnlikelihood for a given correlation matrix

o = ones(samples.npoint,1);
z = zeros(samples.npoint * samples.ndim,1);
one = [o;z];

U = chol(R); % Cholesky Factorisation

mu = (one'*(U\(U'\samples.output_aug)))/(one'*(U\(U'\one)));
sighat=((samples.output_aug-one*mu)'*(U\(U'\(samples.output_aug-one*mu))))/samples.npoint_gek;

% Sum lns of diagonal to find ln(det(R))
LnDetR=2*sum(log(abs(diag(U))));
% ln likelihood is:
lnlikelihood = 0.5*(-samples.npoint_gek*log(sighat)-LnDetR);

end