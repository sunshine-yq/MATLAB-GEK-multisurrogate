function [sample, param] = read_io(options)
% Read the SU2 input and output files 
% These are the GEK samples

%% Supress Warnings
% It's possible that we have some sample points with the same x-y but
% different SA, since we do adjoint on the closest mesh in SU2 some points
% might fall on the same mesh node.
% So suppress warning about duplicate point in XY in interpolation.
warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');

%% Create parameter struct and set integers to each parameter
param.cb1=1; param.sig=2; param.cb2=3; param.kar=4;
param.cw2=5; param.cw3=6; param.cv1=7; param.x=8; param.y=9;

%% Create sample struct to store samples, and specify number of dimensions
sample = struct;
sample.ndim = 9;
sample.nfiles = options.nfiles;

%% Create baseline samples if nfiles specified as 0
if sample.nfiles == 0
    sample = baseline_samples(param, sample, options);
    return
end

%% Read Sample Space as GEK input
% Read samples from files in folder
samfolder = 'Samples/SU2_Input';
allfiles = dir(fullfile(samfolder,'*.dat'));
% Only consider number of files specified
samfiles = allfiles(1:sample.nfiles);
samplearray = cell(length(samfiles),sample.ndim);

formatSpec = '%f%f%f%f%f%f%f%f%f';
for i = 1:length(samfiles)
    fileID = fopen(fullfile(samfiles(i).folder,samfiles(i).name),'r');
    samplearray(i,:)= textscan(fileID, formatSpec,'HeaderLines',1,'Delimiter',',');
    fclose(fileID);
end

% store sample values in variable
sample.input(:,param.cb1) = vertcat(samplearray{:,1});
sample.input(:,param.sig) = vertcat(samplearray{:,2});
sample.input(:,param.cb2) = vertcat(samplearray{:,3});
sample.input(:,param.kar) = vertcat(samplearray{:,4});
sample.input(:,param.cw2) = vertcat(samplearray{:,5});
sample.input(:,param.cw3) = vertcat(samplearray{:,6});
sample.input(:,param.cv1) = vertcat(samplearray{:,7});
sample.input(:,param.x)   = vertcat(samplearray{:,8});
sample.input(:,param.y)   = vertcat(samplearray{:,9});

% Obtain number of Samples
sample.npoint = length(sample.input);

%% Read SU2 results as GEK output to be predicted
% the SU2 output file contains the value of the objective function and the
% value of the gradients, along with the XY coords of the mesh node closest
% to the baseline LHS samples, on which SU2 has evaluated the objective function

% Read output from files in folder
outfolder = 'Samples/SU2_Output';
allfiles = dir(fullfile(outfolder,'*.dat'));
% Only consider number of files specified
outfiles = allfiles(1:sample.nfiles);
% size of SU2 output is 12 columns
outputarray = cell(length(samfiles),sample.ndim+3); 


formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f';
for i = 1:length(outfiles)
    fileID = fopen(fullfile(outfiles(i).folder,outfiles(i).name),'r');
    outputarray(i,:)= textscan(fileID, formatSpec,'HeaderLines',1,'Delimiter','');
    fclose(fileID);
end

% store outputs in one variable
sample.output = vertcat(outputarray{:,3});
sample.outputgrad(:,param.cb1) = vertcat(outputarray{:,4});
sample.outputgrad(:,param.sig) = vertcat(outputarray{:,5});
sample.outputgrad(:,param.cb2) = vertcat(outputarray{:,6});
sample.outputgrad(:,param.kar) = vertcat(outputarray{:,7});
sample.outputgrad(:,param.cw2) = vertcat(outputarray{:,8});
sample.outputgrad(:,param.cw3) = vertcat(outputarray{:,9});
sample.outputgrad(:,param.cv1) = vertcat(outputarray{:,10});
sample.outputgrad(:,param.x)   = vertcat(outputarray{:,11});
sample.outputgrad(:,param.y)   = vertcat(outputarray{:,12});

%% Update input X-Y
% Since the input sample xy location is changed to the nearest mesh node
% inside SU2, update the xy values from the SU2 output file
sample.input(:,param.x) = vertcat(outputarray{:,1});
sample.input(:,param.y) = vertcat(outputarray{:,2});

%% Only retain samples inside global window
inp_aug = [];
out_aug = [];
grd_aug = [];

for i=1:sample.npoint
    if sample.input(i,param.x) >= options.globalx(1) && ...
       sample.input(i,param.x) <= options.globalx(2) && ...
       sample.input(i,param.y) >= options.globaly(1) && ...
       sample.input(i,param.y) <= options.globaly(2)
        inp_aug = cat(1,inp_aug, sample.input(i,:));
        out_aug = cat(1,out_aug, sample.output(i,:));
        grd_aug = cat(1,grd_aug, sample.outputgrad(i,:));
    end
end

sample.input = inp_aug;
sample.output = out_aug;
sample.outputgrad = grd_aug;
sample.npoint = size(inp_aug,1);

%% Augment sample output for GEK
% GEK sample size
sample.npoint_gek = (sample.ndim + 1)*sample.npoint;

% Normalise the output by it's mean and standard deviation
sample.output_mean = mean(sample.output);
sample.output_sd = std(sample.output);
sample.output_norm = (sample.output - sample.output_mean)/sample.output_sd;

% Normalise output gradient by standard deviation of output value
sample.outputgrad_norm = sample.outputgrad/sample.output_sd;

% Find average of abs of gradients to determine most important parameters
% Only done to check later with optimum theta. Not used in GEK
sample.outputgrad_avg = zeros(sample.ndim,1);
for z = 1:sample.ndim
   sample.outputgrad_avg(z) = mean(abs(sample.outputgrad(:,z))); 
end

% Create augmented output matrix by combining output and outputgrad
sample.output_aug = sample.output_norm;
for z = 1:sample.ndim
    sample.output_aug = [sample.output_aug ; sample.outputgrad_norm(:,z)];
end

%% Print header to screen
fprintf('\n***** GEK Problem Definition *****\n');
fprintf('Number of dimensions   = %i\n',sample.ndim);
fprintf('Number of samples      = %i\n', sample.npoint);
fprintf('Number of sample files = %i\n', sample.nfiles);
fprintf('Number of pred points  = %i\n', options.npred);
fprintf('Program objective      = %s\n', options.objective);
if strcmp(options.objective,'batch')
    fprintf('Points in next batch   = %i\n', options.nbatch);
end

end

