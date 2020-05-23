function [] = plotgek(sample, param, pred, batch, pool, options)
% Generate the plots

close all;

% Depending on objective, choose what to plot
if strcmp(options.objective, 'batch')
    plot_mse(sample, param, pred, batch, pool, options);
elseif strcmp(options.objective, 'verify')
    plot_vel(sample, param, pred, options);
end
end

%% MSE Plot
function [] = plot_mse(sample, param, pred, batch, pool, options)
% Plot MSE values of the prediction

% General computations before plotting

% Get the boundaries of the design parameters for plotting
boundary = get_boundary(param, options);
% Get hump surface
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;

% Extract sample points which fall in batch boundaries
sample.input_pool = [];
for i=1:sample.npoint
    if sample.input(i,param.x) >= options.batchxbound(1) && ...
       sample.input(i,param.x) <= options.batchxbound(2) && ...
       sample.input(i,param.y) >= options.batchybound(1) && ...
       sample.input(i,param.y) <= options.batchybound(2)    
       
       sample.input_pool = cat(1,sample.input_pool,sample.input(i,:));
        
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% main figure window
fig = figure;
sgtitle('GEK Prediction Mean Square Error');
addToolbarExplorationButtons(fig);

% MSE of prediction in X-Y space for both pred and pool

% Interpolate the mse
% Pred points
interpx_pred = pred.mapped(:,param.x);
interpy_pred = pred.mapped(:,param.y);
interpz_pred = pred.mse;
interp = scatteredInterpolant(interpx_pred, interpy_pred, interpz_pred, 'linear', 'nearest');
x = linspace(boundary(param.x,1),boundary(param.x,2),1000);
y = linspace(boundary(param.y,1),boundary(param.y,2),1000);
[Xpred,Ypred] = meshgrid(x,y);
% Restack meshgrid to remove y points inside hump
for i=1:length(x)
   if Xpred(1,i) > 0 && Xpred(1,i) < 1
       Ypred(:,i) = linspace(hump_surface(Xpred(1,i)),boundary(param.y,2),length(y))';
   end
end
Zpred = interp(Xpred,Ypred);

% Pool points
interpx_pool = pool.mapped(:,param.x);
interpy_pool = pool.mapped(:,param.y);
interpz_pool = pool.mse;
interp = scatteredInterpolant(interpx_pool, interpy_pool, interpz_pool, 'linear', 'nearest');
x = linspace(options.batchxbound(1),options.batchxbound(2),1000);
y = linspace(options.batchybound(1),options.batchybound(2),1000);
[Xpool,Ypool] = meshgrid(x,y);
% Restack meshgrid to remove y points inside hump 
for i=1:length(x)
   if Xpool(1,i) > 0 && Xpool(1,i) < 1
       if Ypool(1,i) < hump_surface(Xpool(1,i)) % only if hump located in pool window
           Ypool(:,i) = linspace(hump_surface(Xpool(1,i)),options.batchybound(2),length(y))';
       end
   end
end
Zpool = interp(Xpool,Ypool);

% plot the pred mse
p{1} = subplot(3,2,[1,2]);
contlevels = linspace(0,pred.mse_sortval(1),40);
contourf(Xpred,Ypred,Zpred,contlevels,'LineColor','none','HandleVisibility','off');
axis equal; hold on
c{1} = colorbar(p{1});
xlabel('x/c'); ylabel('y/c')
% caxis([0 pred.mse_sortval(1)]);
xlim(boundary(param.x,:));
ylim(boundary(param.y,:));
title('Full X-Y Design Space');

% plot the pool mse
p{2} = subplot(3,2,[3,4]);
contlevels = linspace(0,pool.mse_sortval(1),40);
contourf(Xpool,Ypool,Zpool,contlevels,'LineColor','none','HandleVisibility','off');
axis equal; hold on
c{2} = colorbar(p{2});
xlabel('x/c'); ylabel('y/c')
% caxis([0 pool.mse_sortval(1)]);
xlim(options.batchxbound);
ylim(options.batchybound);
title('Windowed X-Y Design Space');

% plot the batch window rectangle on pred plot
rec_w = options.batchxbound(2) - options.batchxbound(1);
rec_h = options.batchybound(2) - options.batchybound(1);
rectangle(p{1},'Position',[options.batchxbound(1) options.batchybound(1) ...
    rec_w rec_h],'EdgeColor','r','LineWidth',2)

% plot hump and samples and batch on both plots
x = linspace(0,1,1000)';
y = hump_surface(x);
for i=1:length(p)
    area(p{i},x,y,0,'FaceColor','none','HandleVisibility','off')
    plot(p{i},batch.point(:,param.x),batch.point(:,param.y),'*r','linewidth',1)
%     viscircles(p{i},batch.pointxy,batch.radius,'color','k','linewidth',1);
    if i == 1
        plot(p{i},sample.input(:,param.x),sample.input(:,param.y),'xy','linewidth',1);
%         plot(p{i},interpx_pred,interpy_pred,'.m');
    elseif i ==2 % in this plot only plot sample points inside the pool boundary
        plot(p{i},sample.input_pool(:,param.x),sample.input_pool(:,param.y),'xy','linewidth',1);
%         plot(p{i},interpx_pool,interpy_pool,'.m');
    end
end

% Add legend
l = legend(p{1},'Batch','Samples');
l.Position(1) = c{1}.Position(1);
l.Position(2) = c{1}.Position(2) - 0.1;
l.Color = 'k'; l.TextColor = 'w';
l.LineWidth = 1.0; l.FontSize = 9.0; l.FontWeight='bold';

%##########################################################################

% MSE in kar-cb1 space

% Interpolate the mse, full design space
interpx = pred.mapped(:,param.kar);
interpy = pred.mapped(:,param.cb1);
interpz = pred.mse;
interp = scatteredInterpolant(interpx, interpy, interpz, 'linear', 'nearest');

x = linspace(boundary(param.kar,1),boundary(param.kar,2),1000);
y = linspace(boundary(param.cb1,1),boundary(param.cb1,2),1000);
[X,Y] = meshgrid(x,y);
Z = interp(X,Y);

% plot the mse
subplot(3,2,5);
contourf(X,Y,Z,40,'LineColor','none','HandleVisibility','off')
colorbar; hold on
xlabel('kar'); ylabel('cb1')
caxis([0 pred.mse_sortval(1)]);
title('Full kar-cb1 Design Space');

% plot samples and batch
plot(sample.input(:,param.kar),sample.input(:,param.cb1),'xy','linewidth',1);
plot(batch.point(:,param.kar),batch.point(:,param.cb1),'*r','linewidth',1)
% plot(interpx,interpy,'.m');

%##########################################################################

% MSE in sig-cw2 space

% Interpolate the mse, full design space
interpx = pred.mapped(:,param.sig);
interpy = pred.mapped(:,param.cw2);
interpz = pred.mse;
interp = scatteredInterpolant(interpx, interpy, interpz, 'linear', 'nearest');

x = linspace(boundary(param.sig,1),boundary(param.sig,2),1000);
y = linspace(boundary(param.cw2,1),boundary(param.cw2,2),1000);
[X,Y] = meshgrid(x,y);
Z=interp(X,Y);

% plot the mse
subplot(3,2,6);
contourf(X,Y,Z,40,'LineColor','none','HandleVisibility','off')
colorbar; hold on
xlabel('sig'); ylabel('cw2')
caxis([0 pred.mse_sortval(1)]);
title('Full sig-cw2 Design Space');

% plot samples and batch
plot(sample.input(:,param.sig),sample.input(:,param.cw2),'xy','linewidth',1);
plot(batch.point(:,param.sig),batch.point(:,param.cw2),'*r','linewidth',1)
% plot(interpx,interpy,'.m');

%##########################################################################
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