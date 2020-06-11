function [raw_samples, varargout] = revmap_samples(param, mapped_samples, varargin)
% Reverse map samples from design space boundaries back to [0,1]

% number of parameters
ndim = length(fieldnames(param));
% number of samples
nsam = length(mapped_samples);
% Get the boundaries of the parameters
boundary = get_boundary(param);

%% Map the Samples to defined bounds

% Linear mapping of samples to limit bounds
% Map between a & b. a<b

raw_samples = nan(size(mapped_samples));
for j=1:ndim
    a = boundary(j,1);
    b = boundary(j,2);
    for i=1:nsam
        raw_samples(i,j) = (mapped_samples(i,j)-a)/(b-a);
    end
end

%% Map the gradients as well if specified in the input/ouput

if nargin == 3
    varargout = cell(1,1);
    varargout{1} = nan(size(mapped_samples));
    for j=1:ndim
        a = boundary(j,1);
        b = boundary(j,2);
        for i=1:nsam
            varargout{1}(i,j) = varargin{1}(i,j)*(b-a);
        end
    end
end

end

