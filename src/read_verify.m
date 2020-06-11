function [predictions, verify_points] = read_verify(param, options)
% Read the SU2 verify input and output files.
% These are the points against which GEK will be verified

%% Read verify inputs
% Read samples from files in folder
infolder   = sprintf('Samples/Verify');
infile     = sprintf('samples_verify.dat');
fullinfile = fullfile(infolder,infile);

formatSpec = '%f%f%f%f%f%f%f';
fileID = fopen(fullinfile,'r');
samplearray= textscan(fileID, formatSpec,'HeaderLines',1,'Delimiter',',');
fclose(fileID);

% store sample values in variable
verify_points.input(:,param.cb1) = vertcat(samplearray{:,1});
verify_points.input(:,param.sig) = vertcat(samplearray{:,2});
verify_points.input(:,param.cb2) = vertcat(samplearray{:,3});
verify_points.input(:,param.kar) = vertcat(samplearray{:,4});
verify_points.input(:,param.cw2) = vertcat(samplearray{:,5});
verify_points.input(:,param.cw3) = vertcat(samplearray{:,6});
verify_points.input(:,param.cv1) = vertcat(samplearray{:,7});

% Obtain number of Samples
verify_points.npoint = length(verify_points.input);

%% Read verify results
% the SU2 output file contains the value of the objective function 

% Read output from files in folder
outfolder   = sprintf('Samples/Verify');
outfile     = strcat('results_',sprintf('M%.2i',options.activesrrgt),'_verify.dat');
fulloutfile = fullfile(outfolder,outfile);

formatSpec = '%f%f%f';
fileID = fopen(fulloutfile,'r');
outputarray = textscan(fileID, formatSpec,'HeaderLines',1,'Delimiter','');
fclose(fileID);

% store outputs in one variable
verify_points.output = vertcat(outputarray{:,3});

%% Assign predictions points for GEK
% These should be the same as input of the verify points to make comparison
% with verify space consistent
predictions.mapped = verify_points.input;
predictions.npoint = size(predictions.mapped,1);

% Map back to [0,1]
predictions.raw = revmap_samples(param, predictions.mapped);

end

