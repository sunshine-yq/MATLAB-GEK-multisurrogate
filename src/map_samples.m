function [mapped_samples] = map_samples(param,raw_samples, options)
%Function to map samples between [0,1] to boundaries of the
%design parameters. For NASA Hump case.

% number of parameters
ndim = length(fieldnames(param));
% number of samples
nsam = length(raw_samples);
% Get the boundaries of the parameters
boundary = get_boundary(param);

%% Map the Samples to defined bounds

% Linear mapping of samples to limit bounds
% Map between a & b. a<b

mapped_samples = nan(size(raw_samples));
for j=1:ndim
    a = boundary(j,1);
    b = boundary(j,2);
    for i=1:nsam
        mapped_samples(i,j) = raw_samples(i,j)*(b-a)+a;
    end    
end

end

