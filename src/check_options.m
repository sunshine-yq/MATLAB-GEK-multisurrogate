function [options] = check_options(options)
% Assign defaults to empty values if needed and check options entered correctly

% Set defaults if values not chosen by user
if isempty(options.globalx), options.globalx = [-1.0 2.0]; end
if isempty(options.globaly), options.globaly = [ 0.0 0.5]; end
if isempty(options.batchxbound), options.batchxbound = options.globalx; end
if isempty(options.batchybound), options.batchybound = options.globaly; end


% Check if string options entered correctly
if ~strcmp(options.platform,'local') && ~strcmp(options.platform,'iridis')
    error('Wrong platform option specified');
elseif ~strcmp(options.objective,'batch') && ~strcmp(options.objective,'verify')
    error('Wrong objective option specified');
end

% Check if batch boundaries fall outside global boundaries
if options.batchxbound(1)<options.globalx(1) || options.batchxbound(2)>options.globalx(2)
    error('Batch X boundary falls outside of global boundary');
elseif options.batchybound(1)<options.globaly(1) || options.batchybound(2)>options.globaly(2)
    error('Batch X boundary falls outside of global boundary');
end
    
end

