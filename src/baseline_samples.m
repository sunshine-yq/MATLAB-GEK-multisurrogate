function [samples] = baseline_samples(param, options)
% Create baseline samples for SA design space

npoint = options.nsamplesnew;
ndim = numel(fieldnames(param));

skip = floor(rand*1e7);
leap = (nthprime(ndim+1)-1);
halton = haltonset(ndim,'Skip',skip,'Leap',leap);
halton = scramble(halton,'RR2');
raw = net(halton, npoint);

%% Map Samples onto global boundaries
[samples] = map_samples(param, raw, options);

%% Load XY location of surrogate models
surrogates_xy = load('surrogates_xy.mat');
surrogates_xy = surrogates_xy.surrogates_xy;

%% Plots
fig = figure;
sgtitle('Baselines Samples');
addToolbarExplorationButtons(fig);

% sig cw2 samples
subplot(2,1,1)
plot(samples(:,param.sig),samples(:,param.cw2),'*r');
axis equal; grid;
xlabel('sig'); ylabel('cw2');

% kar cb1 samples
subplot(2,1,2)
plot(samples(:,param.kar),samples(:,param.cb1),'*r');
axis equal; grid;
xlabel('kar'); ylabel('cb1');

%% Save samples to csv file for all surrogate models
if options.writetofile
    for ii = 1:options.nmodel
        filename = sprintf('Samples/SU2_Input/I0/samples_M%.2i.dat',ii);
        file = fopen(filename,'w');
        fprintf(file, '%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s \n', ...
            'cb1','sig','cb2','kar','cw2','cw3','cv1','X','Y');
        
        samples_cat = cat(2,samples,repmat(surrogates_xy(ii,:),npoint,1));
        for jj = 1:npoint
            fprintf(file, '%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f \n', ...
                samples_cat(jj,:));
        end
        fclose(file);
    end
end

end

