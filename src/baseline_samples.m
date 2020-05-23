function [sample] = baseline_samples(param, sample, options)
% Create baseline samples for SA & XY design space

sample.npoint = options.nbatch;

skip = floor(rand*1e7);
leap = (nthprime(sample.ndim+1)-1);
halton = haltonset(sample.ndim,'Skip',skip,'Leap',leap);
halton = scramble(halton,'RR2');
sample.raw = net(halton, sample.npoint);

%% Map Samples onto global boundaries
[sample.mapped] = map_samples(param, sample.raw, options);
boundary = get_boundary(param, options);

%% Plots
fig = figure;
sgtitle('Baselines Samples');
addToolbarExplorationButtons(fig);

% XY samples
subplot(2,1,1)
plot(sample.mapped(:,param.x),sample.mapped(:,param.y),'*r');
axis equal;
grid;
hold on;

% Hump surface
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;
x = linspace(0,1,1000)';
y = hump_surface(x);
area(x,y,0,'FaceColor','none','HandleVisibility','off')

xlabel('x/c'); ylabel('y/c');
xlim(boundary(param.x,:));
ylim(boundary(param.y,:));

% kar cb1 samples
subplot(2,1,2)
plot(sample.mapped(:,param.kar),sample.mapped(:,param.cb1),'*r');
axis equal; grid;
xlabel('kar'); ylabel('cb1');

%% Save samples to csv file
if options.writebatch
    file = fopen('Samples/SU2_Input/samples01.dat','w');
    fprintf(file, '%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s \n', ...
        'cb1','sig','cb2','kar','cw2','cw3','cv1','X','Y');
    for i = 1:sample.npoint
        fprintf(file, '%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f \n', ...
            sample.mapped(i,:));
    end
end

end

