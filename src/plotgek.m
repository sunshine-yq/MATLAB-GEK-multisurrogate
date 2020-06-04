function [] = plotgek(samples, param, predictions, nextsamples, options)
% Generate the plots

close all;

% Depending on objective, choose what to plot
if strcmp(options.objective, 'iterate')
    plot_mse(samples, param, predictions, nextsamples, options);
elseif strcmp(options.objective, 'verify')
    plot_vel(samples, param, predictions, options);
end
end

%% MSE Plot
function [] = plot_mse(samples, param, predictions, nextsamples, options)
% Plot MSE values of the prediction. 4 subplots for all 7 SA parameters

% Get the boundaries of the design parameters for plotting
boundary = get_boundary(param);

% Set the pairs to be plotted against eachother in each subplot
plotpairs(1,:) = [param.kar param.cb1];
plotpairs(2,:) = [param.sig param.cw2];
plotpairs(3,:) = [param.cw3 param.cv1];
plotpairs(4,:) = [param.cw2 param.cb2];

% Extract field names of param for axis labels
paramnames = fieldnames(param);

% Plot figures
fig = figure;
sgtitle(sprintf('GEK Prediction MSE - Surrogate M%.2i',options.activesrrgt));
addToolbarExplorationButtons(fig);

for i=1:4    
    % Interpolate the mse
    interpx = predictions.mapped(:,plotpairs(i,1));
    interpy = predictions.mapped(:,plotpairs(i,2));
    interpz = predictions.mse;
    interp = scatteredInterpolant(interpx, interpy, interpz, 'linear', 'nearest');
    
    x = linspace(boundary(plotpairs(i,1),1),boundary(plotpairs(i,1),2),1000);
    y = linspace(boundary(plotpairs(i,2),1),boundary(plotpairs(i,2),2),1000);
    [X,Y] = meshgrid(x,y);
    Z = interp(X,Y);
    
    % plot the mse
    p = subplot(2,2,i);
    contlevels = linspace(0,predictions.mse_sortval(1),40);
    contourf(X,Y,Z,contlevels,'LineColor','none','HandleVisibility','off');
    colorbar; hold on
    xlabel(paramnames(plotpairs(i,1))); ylabel(paramnames(plotpairs(i,2)));
    p.FontWeight = 'bold';
    axis equal
    
    % plot current and new samples
    plot(samples.input(:,plotpairs(i,1)),samples.input(:,plotpairs(i,2)),'xy','linewidth',1);
    plot(nextsamples.input(:,plotpairs(i,1)),nextsamples.input(:,plotpairs(i,2)),'*r','linewidth',1)
%     plot(interpx,interpy,'.m');
end

    l = legend('current','new');
    l.Color = 'k'; l.TextColor = 'w';
    l.LineWidth = 1.0; l.FontSize = 9.0; l.FontWeight='bold';
    l.Position = [0.829 0.927 0.145 0.045];
end

%% Velocity Plot
function [] = plot_vel(sample, param, pred, options)
% Plot value of prediction which is the velocity objective function

% Get the boundaries of the design parameters for plotting
boundary = get_boundary(param, options);
% Get hump surface
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;

% Main figure window
fig = figure(1);
sgtitle('Velocity Objective Function Nominal SA');
addToolbarExplorationButtons(fig);

% output of prediction in X-Y space
p{1} = subplot(3,1,1);

% interpolate the output
interpx = pred.mapped(:,param.x);
interpy = pred.mapped(:,param.y);
interpz = pred.output;
interp = scatteredInterpolant(interpx, interpy, interpz, 'linear', 'nearest');

x = linspace(boundary(param.x,1),boundary(param.x,2),1000);
y = linspace(boundary(param.y,1),boundary(param.y,2),1000);
[X,Y] = meshgrid(x,y);
% Restack meshgrid to remove y points inside hump
for i=1:length(x)
    if X(1,i) > 0 && X(1,i) < 1
        Y(:,i) = linspace(hump_surface(X(1,i)),boundary(param.y,2),length(y))';
    end
end

Z = interp(X,Y);

% plot the output
contourf(X,Y,Z,30,'LineColor','none','HandleVisibility','off')
axis equal; colorbar; hold on; caxis([-1 1.1]);
colormap(p{1},'default');
xlim(boundary(param.x,:));
ylim(boundary(param.y,:));

% plot sample points
plot(sample.input(:,param.x),sample.input(:,param.y),'wx','linewidth',0.8);

title('GEK Prediction');
xlabel('x/c'); ylabel('y/c')

%##########################################################################

% RANS results for Nominal SA
p{2} = subplot(3,1,2);

rans = load('rans.mat');
rans = rans.rans;
rans_velxinterp = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,3), 'linear', 'nearest');
rans_velyinterp = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,4), 'linear', 'nearest');

% Acquire velocities from interpolated function
rans_velx = rans_velxinterp(X,Y);
rans_vely = rans_velyinterp(X,Y);

% calculate objective function
rans_velmag = sqrt(rans_velx.^2 + rans_vely.^2);
rans_velang = atan2(rans_vely, rans_velx);
rans_obj = rans_velmag .* rans_velang;

contourf(X,Y,rans_obj,30,'LineColor','none','HandleVisibility','off')
axis equal; colorbar; hold on; caxis([-1 1.1]);
colormap(p{2},'default');
xlim(boundary(param.x,:));
ylim(boundary(param.y,:));

title('RANS SU2');
xlabel('x/c'); ylabel('y/c')

%##########################################################################

% Difference between RANS and GEK
p{3} = subplot(3,1,3);

diff = abs(Z - rans_obj);

contlevels = linspace(0,1,20);
contourf(X,Y,diff,contlevels,'LineColor','none','HandleVisibility','off')
axis equal; colorbar; hold on;
colormap(p{3},'jet');
xlim(boundary(param.x,:));
ylim(boundary(param.y,:));

% plot sample points
plot(sample.input(:,param.x),sample.input(:,param.y),'wx','linewidth',0.8);

title('|RANS - GEK|');
xlabel('x/c'); ylabel('y/c')
xlim(boundary(param.x,:));
ylim(boundary(param.y,:));

%##########################################################################

% Figure 2, plot MSE of the prediction with Nominal SA

% Use previous meshgrid and interpx and interpy
interpz = pred.mse;
interp = scatteredInterpolant(interpx, interpy, interpz, 'linear', 'nearest');
Z = interp(X,Y);

% Main figure window
fig = figure(2);
addToolbarExplorationButtons(fig);

% plot the MSE
contourf(X,Y,Z,30,'LineColor','none','HandleVisibility','off')
axis equal; colorbar; hold on;
caxis([0 pred.mse_sortval(1)]);
xlim(boundary(param.x,:));
ylim(boundary(param.y,:));

% plot sample points
plot(sample.input(:,param.x),sample.input(:,param.y),'rx','linewidth',0.8);

xlabel('x/c'); ylabel('y/c')
title('GEK Prediction MSE Nominal SA');
p{4} = fig.CurrentAxes;

%##########################################################################

% plot hump
x = linspace(0,1,1000)';
y = hump_surface(x);
for i = 1:length(p)
    area(p{i},x,y,0,'FaceColor','none','HandleVisibility','off')
end

%##########################################################################
end