function [boundary] = get_boundary(param)
% Define the physical boundary of the design parameters

boundary = zeros(numel(fieldnames(param)),2);
boundary(param.cb1,:) = [0.129 0.14];
boundary(param.sig,:) = [0.6 1.4];
boundary(param.cb2,:) = [0.61 0.7];
boundary(param.kar,:) = [0.36 0.42];
boundary(param.cw2,:) = [0.055 0.353];
boundary(param.cw3,:) = [1.5 2.75];
boundary(param.cv1,:) = [6.9 7.5];

end

