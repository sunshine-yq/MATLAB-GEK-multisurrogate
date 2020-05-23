%% Code to compare the velocity field between nominal SA SU2 and LES and Experiment
% Read results from the nominal steady SU2 run and compare
% Since meshes are different we need to interpolate

clear; close all
warning('off','all')

%% Read results and do interpolation
les = load('les.mat'); les = les.les;
rans = load('rans.mat'); rans = rans.rans;

% Do interpolation since meshes are different and we can't compare node to node
rans_velxint = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,3), 'linear', 'nearest');
rans_velyint = scatteredInterpolant(rans(:,1),rans(:,2),rans(:,4), 'linear', 'nearest');
les_velxint = scatteredInterpolant(les(:,1),les(:,2),les(:,3), 'linear', 'nearest');
les_velyint = scatteredInterpolant(les(:,1),les(:,2),les(:,4), 'linear', 'nearest');

%% Create meshgrid for plotting contours
% Get hump surface
hump_surface = load('hump_surface.mat');
hump_surface = hump_surface.hump_surface;

% Number of points to plot contour
xpoint = 6000;
ypoint = 1000;

% Global limits
xbound = [0.7, 1.5];
ybound = [0, 0.1];
% xbound = [-1, max(les(:,1))];
% ybound = [0, max(les(:,2))];

x = linspace(xbound(1),xbound(2),xpoint)';
y = linspace(ybound(1),ybound(2),ypoint)';
[X,Y] = meshgrid(x,y);

% Restack meshgrid to remove y points inside hump
for i=1:xpoint
   if X(1,i) > 0 && X(1,i) < 1
       Y(:,i) = linspace(hump_surface(X(1,i)),ybound(2),ypoint)';
   end
end

%% Acquire velocities from interpolated function
rans_velx = rans_velxint(X,Y);
rans_vely = rans_velyint(X,Y);
les_velx  = les_velxint(X,Y);
les_vely  = les_velyint(X,Y);

% Find magnitude and angle of velocity vector
rans_velmag = sqrt(rans_velx.^2 + rans_vely.^2);
rans_velang = atan2(rans_vely, rans_velx);
les_velmag  = sqrt(les_velx.^2 + les_vely.^2);
les_velang  = atan2(les_vely, les_velx);

% Find objective function value
rans_objfunc = rans_velmag.*rans_velang;
les_objfunc  = les_velmag.*les_velang;

%% Find difference
diff_velmag  = abs(rans_velmag-les_velmag);
diff_velang  = abs(rans_velang-les_velang);
diff_objfunc = abs(rans_objfunc-les_objfunc);

%% Find locations of surrogate models
nmod = 50; % number of models
maxrad = 0.01; % max radius
ymin = 0.00015; % smallest y allowed (expect highest MSE point)

diff_objfunc_col = reshape(diff_objfunc,[xpoint*ypoint,1]);
xmeshgrid_col = reshape(X,[xpoint*ypoint,1]);
ymeshgrid_col = reshape(Y,[xpoint*ypoint,1]);

[diff_objvalsort, diff_objindsort] = sort(diff_objfunc_col,'descend');

model_coor = [xmeshgrid_col(diff_objindsort(1)), ymeshgrid_col(diff_objindsort(1))];

ii = 2;

while size(model_coor,1) < nmod && ii <= xpoint*ypoint
    candidate = [xmeshgrid_col(diff_objindsort(ii)), ymeshgrid_col(diff_objindsort(ii))];
    qualified = true;
    
    for j=1:size(model_coor,1)
       dist =  norm(candidate - model_coor(j,:));
       
       if dist <= maxrad || candidate(2)<=ymin
           qualified = false;
           break;
       end
    end
    
    if qualified
        model_coor = cat(1,model_coor,candidate);
    end
    
    ii = ii+1;    
end

%% Plot of difference
figure(1);

% difference threshold
thresh = 0;

p{1}=subplot(3,1,1);
levels = linspace(thresh,max(diff_velmag,[],'all'),40);
contourf(X,Y,diff_velmag,levels,'LineColor','none')
axis equal; hold on;
title('Difference in Velocity Magnitude - RANS & LES')
colorbar
colormap('jet')

p{2}=subplot(3,1,2);
levels = linspace(thresh,max(diff_velang,[],'all'),40);
contourf(X,Y,diff_velang,levels,'LineColor','none')
axis equal; hold on;
title('Difference in Velocity Angle - RANS & LES')
colorbar
colormap('jet')

p{3}=subplot(3,1,3);
levels = linspace(thresh,max(diff_objfunc,[],'all'),40);
contourf(X,Y,diff_objfunc,levels,'LineColor','none')
axis equal; hold on;
title('Difference in Velocity ObjFunc - RANS & LES')
colorbar
colormap('jet')

fig = figure(2);
contourf(X,Y,rans_objfunc,40,'LineColor','none')
axis equal; hold on;
title('RANS objfunc')
colorbar
colormap('jet')
p{4} = fig.CurrentAxes;

%% Plot hump and surrogate locations on all figures and set axis labels

xhump = linspace(xbound(1),xbound(2),1000)';
yhump = hump_surface(xhump);

for i=1:length(p)
    area(p{i},xhump,yhump,0,'FaceColor','none')
    plot(p{i},model_coor(:,1),model_coor(:,2),'rx','linewidth',3);
    xlabel('x/c');
    ylabel('y/c');
end
