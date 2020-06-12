function [] = plotgek(samples, param, predictions, nextsamples, verifypoints, options)
% Generate the plots

close all;

% Depending on objective, choose what to plot
if strcmp(options.objective, 'iterate')
    plot_mse(samples, param, predictions, nextsamples, options);
elseif strcmp(options.objective, 'verify')
    plot_ver(samples, param, predictions, verifypoints, options);
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
sgtitle(sprintf('GEK Prediction MSE - Surrogate M%.2i Iteration I%.2i',options.activesrrgt,options.nfiles));
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
    plot(samples.inputmapped(:,plotpairs(i,1)),samples.inputmapped(:,plotpairs(i,2)),'xy','linewidth',1);
    plot(nextsamples.inputmapped(:,plotpairs(i,1)),nextsamples.inputmapped(:,plotpairs(i,2)),'*r','linewidth',1)
    %     plot(interpx,interpy,'.m');
    
end

l = legend('current','new');
l.Color = 'k'; l.TextColor = 'w';
l.LineWidth = 1.0; l.FontSize = 9.0; l.FontWeight='bold';
l.Position = [0.829 0.927 0.145 0.045];
end

%% Verify Plot
function [] = plot_ver(samples, param, predictions, verifypoints, options)
% Plot velocity objective function: prediction vs full order
% Difference contour between GEK and Full Order
% Actual Values as pointcloud

% Get the boundaries of the design parameters for plotting
boundary = get_boundary(param);

% Set the pairs to be plotted against eachother in each subplot
plotpairs(1,:) = [param.kar param.cb1];
plotpairs(2,:) = [param.sig param.cw2];
plotpairs(3,:) = [param.cw3 param.cv1];
plotpairs(4,:) = [param.cw2 param.cb2];

% Extract field names of param for axis labels
paramnames = fieldnames(param);

% Figure 1 difference contour
fig1 = figure(1);
sgtitle(sprintf('Velocity Function: |GEK Prediction - Full Order| - Surrogate M%.2i Iteration I%.2i',options.activesrrgt,options.nfiles));
addToolbarExplorationButtons(fig1);

% Figure 2 actual value pointcloud
fig2 = figure(2);
sgtitle(sprintf('Velocity Function: Absolute Values - Surrogate M%.2i Iteration I%.2i',options.activesrrgt,options.nfiles));
addToolbarExplorationButtons(fig2);

for i=1:4
    % Interpolate the prediction output and full order output
    interpx_gek = predictions.mapped(:,plotpairs(i,1));
    interpy_gek = predictions.mapped(:,plotpairs(i,2));
    interpz_gek = predictions.output;
    interp_gek = scatteredInterpolant(interpx_gek, interpy_gek, interpz_gek, 'natural', 'nearest');
    
    interpx_full = verifypoints.input(:,plotpairs(i,1));
    interpy_full = verifypoints.input(:,plotpairs(i,2));
    interpz_full = verifypoints.output;
    interp_full = scatteredInterpolant(interpx_full, interpy_full, interpz_full, 'natural', 'nearest');
    
    x = linspace(boundary(plotpairs(i,1),1),boundary(plotpairs(i,1),2),1000);
    y = linspace(boundary(plotpairs(i,2),1),boundary(plotpairs(i,2),2),1000);
    [X,Y] = meshgrid(x,y);
    
    Z_gek  = interp_gek(X,Y);
    Z_full = interp_full(X,Y);
    Z = abs(Z_gek-Z_full);
    
    %~~~ FIGURE 1
    % plot the difference
    set(0, 'CurrentFigure', fig1);
    p = subplot(2,2,i);
    contourf(X,Y,Z,20,'LineColor','none','HandleVisibility','off');
    colorbar; hold on
    xlabel(paramnames(plotpairs(i,1))); ylabel(paramnames(plotpairs(i,2)));
    p.FontWeight = 'bold';
    axis equal
    % plot samples
    plot(samples.inputmapped(:,plotpairs(i,1)),samples.inputmapped(:,plotpairs(i,2)),'xr','linewidth',1);
    %     plot(interpx,interpy,'.m');
    
    %~~~ FIGURE 2
    % plot the prediction pointcloud
    set(0, 'CurrentFigure', fig2);
    p = subplot(2,2,i);
    %     mesh(X,Y,Z_gek,'HandleVisibility','off');
    plot3(predictions.mapped(:,plotpairs(i,1)),predictions.mapped(:,plotpairs(i,2)), ...
        predictions.output,'ok','linewidth',2,'MarkerSize',5)
    xlabel(paramnames(plotpairs(i,1))); ylabel(paramnames(plotpairs(i,2)));
    zlabel('Velocity Function');
    p.FontWeight = 'bold';
    title('Prediction');
    hold on; grid
    plot3(verifypoints.input(:,plotpairs(i,1)),verifypoints.input(:,plotpairs(i,2)), ...
        verifypoints.output,'.r','linewidth',2,'MarkerSize',12)
    
end
l = legend('Prediction','Full Order');
l.LineWidth = 1.0; l.FontSize = 9.0; l.FontWeight='bold';
l.Position = [0.829 0.927 0.145 0.045];

end