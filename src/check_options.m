function check_options(options)
% Check options entered correctly

% Check if string options entered correctly
if ~strcmp(options.platform,'local') && ~strcmp(options.platform,'iridis')
    error('Wrong platform option specified');
elseif ~strcmp(options.objective,'iterate') && ~strcmp(options.objective,'verify')
    error('Wrong objective option specified');
end
   
end

