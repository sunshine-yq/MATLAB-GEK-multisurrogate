function [mapped_samples] = map_samples(param,raw_samples, options)
%MAPTOBOUNDS Function to map samples between [0,1] to boundaries of the
%design parameters. For NASA Hump case.

% Input: par = struct with parameter names and assinged integer vals
%        raw_samples = samples between [0,1] obtained using quasi-random sets
% varargin: if required to change the xy mapping boundaries, add them here

% Output: mapped_samples = samples mapped linearly to the defined bounds

% number of parameters
npar = length(fieldnames(param));
% number of samples
nsam = length(raw_samples);
% hump surface
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;

% Get the boundaries of the parameters
boundary = get_boundary(param, options);

%% Map the Samples to defined bounds

% Linear mapping of samples to limit bounds
% Map between a & b. a<b

mapped_samples = nan(size(raw_samples));
for j=1:npar
    if j ~= param.y % not y, map based on fixed bounds
        a = boundary(j,1);
        b = boundary(j,2);
        for i=1:nsam
            mapped_samples(i,j) = raw_samples(i,j)*(b-a)+a;
        end
    else % y: lower bound changes based on hump coordinates
        x = mapped_samples(:,param.x);
        % lower bound "a" should be whichever is the highest b/w hump
        % coordinate and specified lower limit of y boundary. This is
        % important when we have a batch window with it's bottom frame above
        % the hump.
        a = max(hump_surface(x),boundary(param.y,1));
        b = boundary(param.y,2);
        for i=1:nsam       
            mapped_samples(i,j) = raw_samples(i,j)*(b-a(i))+a(i);
        end
    end
end

end

