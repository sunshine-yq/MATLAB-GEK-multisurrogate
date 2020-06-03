function [samples] = read_io(param, options)
% Read the SU2 input and output files. These are the GEK samples

%% Model dimension
samples.ndim = length(fieldnames(param));
%% Read Sample Space as GEK input
% Read samples from files in folder
infolder = strcat('Samples/',sprintf('M%.2i',options.activesrrgt),'/SU2_Input');
allfiles = dir(fullfile(infolder,'*.dat'));
% Only consider number of files specified
samfiles = allfiles(1:options.nfiles);
% SU2 input file has 9 columns
samplearray = cell(length(samfiles),9);

formatSpec = '%f%f%f%f%f%f%f%f%f';
for i = 1:length(samfiles)
    fileID = fopen(fullfile(samfiles(i).folder,samfiles(i).name),'r');
    samplearray(i,:)= textscan(fileID, formatSpec,'HeaderLines',1,'Delimiter',',');
    fclose(fileID);
end

% store sample values in variable
samples.input(:,param.cb1) = vertcat(samplearray{:,1});
samples.input(:,param.sig) = vertcat(samplearray{:,2});
samples.input(:,param.cb2) = vertcat(samplearray{:,3});
samples.input(:,param.kar) = vertcat(samplearray{:,4});
samples.input(:,param.cw2) = vertcat(samplearray{:,5});
samples.input(:,param.cw3) = vertcat(samplearray{:,6});
samples.input(:,param.cv1) = vertcat(samplearray{:,7});

% Obtain number of Samples
samples.npoint = length(samples.input);

%% Read SU2 results as GEK output to be predicted
% the SU2 output file contains the value of the objective function and the
% value of the gradients

% Read output from files in folder
outfolder = strcat('Samples/',sprintf('M%.2i',options.activesrrgt),'/SU2_Output');
allfiles = dir(fullfile(outfolder,'*.dat'));
% Only consider number of files specified
outfiles = allfiles(1:options.nfiles);
% SU2 output file has 12 columns
outputarray = cell(length(samfiles),12); 


formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f';
for i = 1:length(outfiles)
    fileID = fopen(fullfile(outfiles(i).folder,outfiles(i).name),'r');
    outputarray(i,:)= textscan(fileID, formatSpec,'HeaderLines',1,'Delimiter','');
    fclose(fileID);
end

% store outputs in one variable
samples.output                  = vertcat(outputarray{:,3});
samples.outputgrad(:,param.cb1) = vertcat(outputarray{:,4});
samples.outputgrad(:,param.sig) = vertcat(outputarray{:,5});
samples.outputgrad(:,param.cb2) = vertcat(outputarray{:,6});
samples.outputgrad(:,param.kar) = vertcat(outputarray{:,7});
samples.outputgrad(:,param.cw2) = vertcat(outputarray{:,8});
samples.outputgrad(:,param.cw3) = vertcat(outputarray{:,9});
samples.outputgrad(:,param.cv1) = vertcat(outputarray{:,10});

%% Augment sample output for GEK
% GEK sample size
samples.npoint_gek = (samples.ndim + 1)*samples.npoint;

% Normalise the output by it's mean and standard deviation
samples.output_mean = mean(samples.output);
samples.output_sd = std(samples.output);
samples.output_norm = (samples.output - samples.output_mean)/samples.output_sd;

% Normalise output gradient by standard deviation of output value
samples.outputgrad_norm = samples.outputgrad/samples.output_sd;

% Find average of abs of gradients to determine most important parameters
% Only done to check later with optimum theta. Not used in GEK
samples.outputgrad_avg = zeros(samples.ndim,1);
for z = 1:samples.ndim
   samples.outputgrad_avg(z) = mean(abs(samples.outputgrad(:,z))); 
end

% Create augmented output matrix by combining output and outputgrad
samples.output_aug = samples.output_norm;
for z = 1:samples.ndim
    samples.output_aug = [samples.output_aug ; samples.outputgrad_norm(:,z)];
end

%% Print header to screen
% % fprintf('\n***** GEK Problem Definition *****\n');
% % fprintf('Number of dimensions   = %i\n',samples.ndim);
% % fprintf('Number of samples      = %i\n', sample.npoint);
% % fprintf('Number of sample files = %i\n', options.nfiles);
% % fprintf('Number of pred points  = %i\n', options.npred);
% % fprintf('Program objective      = %s\n', options.objective);
% % if strcmp(options.objective,'batch')
% %     fprintf('Points in next batch   = %i\n', options.nbatch);
% % end

end

