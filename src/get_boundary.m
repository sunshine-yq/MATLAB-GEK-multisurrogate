function [boundary] = get_boundary(par, options)
% define the physical boundary of the design parameters
% input: struct par containing the parameters, options
% output: boundaries for each parameter

boundary = zeros(numel(fieldnames(par)),2);
boundary(par.cb1,:) = [0.129 0.14];
boundary(par.sig,:) = [0.6 1.4];
boundary(par.cb2,:) = [0.61 0.7];
boundary(par.kar,:) = [0.36 0.42];
boundary(par.cw2,:) = [0.055 0.353];
boundary(par.cw3,:) = [1.5 2.75];
boundary(par.cv1,:) = [6.9 7.5];
boundary(par.x,:)   = options.globalx;
boundary(par.y,:)   = options.globaly;

end

