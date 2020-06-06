function baseline_samples(param, options)
% Create baseline samples for SA design space

npoint = options.nnextsamples;
ndim   = length(fieldnames(param));
% Halton sequence
skip = floor(rand*1e7);
leap = (nthprime(ndim+1)-1);
halton = haltonset(ndim,'Skip',skip,'Leap',leap);
halton = scramble(halton,'RR2');
rawsamples = net(halton, npoint);

%% Map Samples onto global boundaries
[samples] = map_samples(param, rawsamples);

%% Load XY location of surrogate models
% xy still needs to be specified in SU2 for fixed location of surrogates
surrogates_xy = load('surrogates_xy.mat');
surrogates_xy = surrogates_xy.surrogates_xy;

%% Save samples to csv files for all surrogate models
if options.writetofile
    for ii = 1:options.nsurrogates
        filename = sprintf('Samples/M%.2i/SU2_Input/samples_M%.2i_I01.dat',ii,ii);
        file = fopen(filename,'w');
        fprintf(file, '%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s \n', ...
            'cb1','sig','cb2','kar','cw2','cw3','cv1','X','Y');
        
        samples_xycat = cat(2,samples,repmat(surrogates_xy(ii,:),npoint,1));
        for jj = 1:npoint
            fprintf(file, '%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f \n', ...
                samples_xycat(jj,:));
        end
        fclose(file);
    end
end

%% Plots
fig = figure;
sgtitle('Baselines Samples');
addToolbarExplorationButtons(fig);

% kar cb1 samples
subplot(2,2,1)
plot(samples(:,param.kar),samples(:,param.cb1),'*r');
axis equal; grid;
xlabel('kar'); ylabel('cb1');

% sig cw2 samples
subplot(2,2,2)
plot(samples(:,param.sig),samples(:,param.cw2),'*r');
axis equal; grid;
xlabel('sig'); ylabel('cw2');

% cw3 cv1 samples
subplot(2,2,3)
plot(samples(:,param.cw3),samples(:,param.cv1),'*r');
axis equal; grid;
xlabel('cw3'); ylabel('cv1');

% cw2 cb2 samples
subplot(2,2,4)
plot(samples(:,param.cw2),samples(:,param.cb2),'*r');
axis equal; grid;
xlabel('cw2'); ylabel('cb2');

end

