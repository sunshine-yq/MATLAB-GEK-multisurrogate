function [nextsamples] = next_iteration(predictions, options)
% Determine prediction points with the highest MSE. These become the samples
% for the next SU2 iteration.

% Store highest MSE prediction points as new samples
% only storing mapped samples to pass on to SU2 later
% the indecies of predictions mapped and raw are consistent anyway
nextsamples.inputmapped = predictions.mapped(predictions.mse_sortindex(1:options.nnextsamples),:);
nextsamples.mse   = predictions.mse_sortval(1:options.nnextsamples);

%% Write new samples to file
if options.writetofile
    % Load XY location of surrogate models
    % xy still needs to be specified in SU2 for fixed location of surrogates
    surrogates_xy = load('surrogates_xy.mat');
    surrogates_xy = surrogates_xy.surrogates_xy;
    
    
    nextiterno = options.nfiles + 1; % one after current file number
    folder = strcat('Samples/',sprintf('M%.2i',options.activesrrgt),'/SU2_Input');
    filename = sprintf('samples_M%.2i_I%.2i.dat', options.activesrrgt, nextiterno);
    file = fopen(fullfile(folder,filename),'w');
    fprintf(file, '%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s \n', ...
        'cb1','sig','cb2','kar','cw2','cw3','cv1','X','Y');
    
    samples_xycat = ...
        cat(2,nextsamples.inputmapped,repmat(surrogates_xy(options.activesrrgt,:),options.nnextsamples,1));
    
    for i = 1:options.nnextsamples
        fprintf(file, '%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f \n', ...
            samples_xycat(i,:));
    end
end

end

