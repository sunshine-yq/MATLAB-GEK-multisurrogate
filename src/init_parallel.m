function [] = init_parallel(options)
% Set number of nodes to initialise parallel run

if strcmp(options.platform,'local')
    numcpu = 4;
elseif strcmp(options.platform,'iridis')
    numcpu = 40;
else
    error('Invalid platform name');
end

p = gcp('nocreate');
if isempty(p)
    parpool(numcpu);
end

end
