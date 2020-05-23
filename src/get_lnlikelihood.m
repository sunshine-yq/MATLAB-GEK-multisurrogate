function [lnlikelihood] = get_lnlikelihood(R, sample)
% Function to get the lnlikelihood for a given correlation matrix

o = ones(sample.npoint,1);
z = zeros(sample.npoint * sample.ndim,1);
one = [o;z];

U = chol(R); % Cholesky Factorisation

mu = (one'*(U\(U'\sample.output_aug)))/(one'*(U\(U'\one)));
sighat=((sample.output_aug-one*mu)'*(U\(U'\(sample.output_aug-one*mu))))/sample.npoint_gek;

% Sum lns of diagonal to find ln(det(R))
LnDetR=2*sum(log(abs(diag(U))));
% ln likelihood is:
lnlikelihood = 0.5*(-sample.npoint_gek*log(sighat)-LnDetR);

end