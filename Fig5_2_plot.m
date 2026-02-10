clear;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'data');

dataFile = fullfile(dataDir,'data_Fig5.mat');

if ~exist(dataFile,'file')
    error(['Required data file not found: %s\n' ...
           'Please run Fig5_1_run.m first, or download the precomputed data.'], ...
           dataFile);
end
%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%% 
load(dataFile)
Xv = data_Fig5.grid.X;
Yv = data_Fig5.grid.Y;

%% plot Fig5a-c
figure; 
subplot(1,3,1); 
axis equal;
axis([0 1 0 1]);
axis off;
hold on;
pcolor(Xv, Yv, data_Fig5.Md(:,:,1));shading flat;
colorbar;
title('Fig5a');
subplot(1,3,2); 
axis equal;
axis([0 1 0 1]);
axis off;
hold on;
pcolor(Xv, Yv, data_Fig5.Ms(:,:,2));shading flat;
colorbar;
title('Fig5b');
ax3 = subplot(1,3,3); 
axis equal;
axis([0 1 0 1]);
axis off;
hold on;
dM = data_Fig5.Md(:,:,1)-data_Fig5.Ms(:,:,2);
pcolor(Xv, Yv, dM);shading flat;
colorbar;
MIN_dM = min(dM,[],'all');
MAX_dM = max(dM,[],'all');
CustomColormap = colormapGeneration(MIN_dM, MAX_dM);
colormap(ax3, CustomColormap);
title('Fig5c');

%% plot Fig5d
[xGrid, yGrid] = meshgrid(Xv, Yv);
Tri = delaunay(xGrid, yGrid);
figure; hold on;

% using two colors for upper/lower surfaces
ColorSurf = [1 0.5 0.5; ...
             0.1 0.1 0.9; ...
             0.4 0.4 0.9; ...
             0.7 0.7 0.9];
ColorSurfBack = 1-(1-ColorSurf)*0.5;
dd = 0.001;

% upper surface
trisurf(Tri, xGrid, yGrid, data_Fig5.Md(:,:,1)+dd, 'EdgeColor','None', 'FaceColor',ColorSurf(1,:),'FaceAlpha',1);
trisurf(Tri, xGrid, yGrid, data_Fig5.Ms(:,:,1)+dd, 'EdgeColor','None', 'FaceColor',ColorSurf(2,:),'FaceAlpha',1);
trisurf(Tri, xGrid, yGrid, data_Fig5.Ms(:,:,2)+dd, 'EdgeColor','None', 'FaceColor',ColorSurf(3,:),'FaceAlpha',1);
trisurf(Tri, xGrid, yGrid, data_Fig5.Ms(:,:,3)+dd, 'EdgeColor','None', 'FaceColor',ColorSurf(4,:),'FaceAlpha',1);

% lower surface
trisurf(Tri, xGrid, yGrid, data_Fig5.Md(:,:,1), 'EdgeColor','None', 'FaceColor',ColorSurfBack(1,:),'FaceAlpha',1);
trisurf(Tri, xGrid, yGrid, data_Fig5.Ms(:,:,1), 'EdgeColor','None', 'FaceColor',ColorSurfBack(2,:),'FaceAlpha',1);
trisurf(Tri, xGrid, yGrid, data_Fig5.Ms(:,:,2), 'EdgeColor','None', 'FaceColor',ColorSurfBack(3,:),'FaceAlpha',1);
trisurf(Tri, xGrid, yGrid, data_Fig5.Ms(:,:,3), 'EdgeColor','None', 'FaceColor',ColorSurfBack(4,:),'FaceAlpha',1);

% A point for example
ix = 400; iy = 125;
plot3(Xv(ix), Yv(iy), data_Fig5.Ms(iy,ix,1),'ko','LineWidth',1,'MarkerFaceColor','k','MarkerSize',3);
plot3(Xv(ix), Yv(iy), data_Fig5.Ms(iy,ix,2),'ko','LineWidth',1,'MarkerFaceColor','k','MarkerSize',3);
plot3(Xv(ix), Yv(iy), data_Fig5.Ms(iy,ix,3),'ko','LineWidth',1,'MarkerFaceColor','k','MarkerSize',3);
plot3(Xv(ix), Yv(iy), data_Fig5.Md(iy,ix,1),'ko','LineWidth',1,'MarkerFaceColor','k','MarkerSize',3);
plot3(Xv(ix), Yv(iy), 0,'ko','LineWidth',1,'MarkerFaceColor','k','MarkerSize',6);

% Distances to three edge
d1 = abs(sqrt(3)*Xv(ix)/2-1/2*Yv(iy));
d2 = abs(sqrt(3)*Xv(ix)/2+1/2*Yv(iy)-sqrt(3)/2);
plot([Xv(ix) Xv(ix)], [0 Yv(iy)],'k--','LineWidth',1);
plot([Xv(ix) Xv(ix)-sqrt(3)/2*d1], [Yv(iy) Yv(iy)+d1/2],'k--','LineWidth',1);
plot([Xv(ix) Xv(ix)+sqrt(3)/2*d2], [Yv(iy) Yv(iy)+d2/2],'k--','LineWidth',1);

L = 1000;
% Inside boundary of FDOL (piecewise planar structure for parameters used in Fig5d)
x1 = linspace(2/9,11/18,L); y1 = sqrt(3)*(x1-2/9); 
z1 = max(0,2*(x1-2/9)*3-2/3)+max(0,2*(11/18-x1)*3-2/3);
plot3(x1,y1,z1,'k','LineWidth',1);

x2 = linspace(1/9,8/9,L); y2 = sqrt(3)/9*ones(1,L); 
z2 = max(0,(x2-1/9)*3-2/3)+max(0,(8/9-x2)*3-2/3);
plot3(x2,y2,z2,'k','LineWidth',1);

x3 = linspace(7/18,7/9,L); y3 = 7/18*sqrt(3)-sqrt(3)*(x3-7/18); 
z3 = max(0,2*(x3-7/18)*3-2/3)+max(0,2*(7/9-x3)*3-2/3);
plot3(x3,y3,z3,'k','LineWidth',1);


% Outside boundary
boundaryX = [linspace(0,1/2,L);...
            linspace(0,1,L);...
            linspace(1/2,1,L)];
boundaryY = [sqrt(3)*linspace(0,1/2,L);...
            zeros(1,L);...
            sqrt(3)-sqrt(3)*linspace(1/2,1,L)];

% reference axis
plot3([1,1],[0,0], [0,3],'k','LineWidth',2);
plot3(boundaryX', boundaryY', zeros(size(boundaryX))','k', 'LineWidth',1.5);

title('Fig5d');
axis([0 1 0 1 0 3]);
view([190 14]);

%% Plot Fig5f
gp_5f = 0.4;
ep_5f = 0.4;

n_scan = 1:0.1:10;
g_n = metabolic_cost(n_scan, gp_5f);        % g^[n] = 1 + (n-1) g_p
sigma_n = uptake_efficiency(n_scan, ep_5f); % sigma^[n] = 1 / (1 + e_p (n-1))

figure; hold on;

yyaxis left
plot(n_scan, g_n);
ylim([0 5])
ylabel('Relative metabolic cost g');
xlabel('Number of pathways n');

yyaxis right
plot(n_scan,sigma_n);
ylim([0 1])
ylabel('Relative pathway efficiency');
box on;
title('Fig5f');

%% Plot Fig5g
gp_scan = 0:0.01:2;
ep_boundary = (1-gp_scan)./(1+gp_scan);

figure; hold on;
plot(gp_scan, ep_boundary, 'k');
axis([0 2 0 2]);
xlabel('g_p');
ylabel('e_p');
box on;
title('Fig5g');
%% Plot Fig5h
NT_5h = 4;
ep_5h = 0.1; gp_5h = 0.1;

% Ideal curve for n (not trunctated by NT)
mu0_min = mu0_min_optimum(gp_5h, ep_5h);
mu0_scan = exp(log(mu0_min): 0.01: log(100));

np = max(np_continuous(mu0_scan, gp_5h, ep_5h), 1);

% The curve in practice (trunctated by NT)
mu0_min2 = mu0_min_SS(NT_5h, gp_5h, ep_5h);
mu0_scan2 = exp(log(mu0_min2): 0.01: log(100));
np2 = max(np_continuous(mu0_scan2, gp_5h, ep_5h), 1);
np2(np2>NT_5h) = NT_5h;

figure; 
semilogx(mu0_scan, np); hold on;
semilogx(mu0_scan2, np2)
xlim([0.1 100])
ylim([0 10])
xticks([]);
yticks([]);
box on;
title('Fig5h');

%% Plot Fig5i
NT_5i = 10;
ep_5i = 0.1; gp_5i = 0.1;

% Ideal curve for n (not trunctated by NT)
mu0_min = mu0_min_optimum(gp_5i, ep_5i);
mu0_scan = exp(log(mu0_min): 0.01: log(100));

np = max(np_continuous(mu0_scan, gp_5i, ep_5i), 1);

% The curve in practice (trunctated by NT)
mu0_min2 = mu0_min_SS(NT_5i, gp_5i, ep_5i);
mu0_scan2 = exp(log(mu0_min2): 0.01: log(100));
np2 = max(np_continuous(mu0_scan2, gp_5i, ep_5i), 1);
np2(np2>NT_5i) = NT_5i;

figure; 
semilogx(mu0_scan, np); hold on;
semilogx(mu0_scan2, np2)
xlim([0.1 100])
ylim([0 10])
xticks([]);
yticks([]);
box on;
title('Fig5i');

%% ===================== Local helper functions ===========================
function g_n = metabolic_cost(n, gp)
    % g^[n] = 1 + (n-1) g_p
    g_n = 1 + gp.*(n - 1);
end

function sigma_n = uptake_efficiency(n, ep)
    % sigma^[n] = 1 / (1 + e_p (n-1))
    sigma_n = 1 ./ (1 + ep.*(n - 1));
end

function np = np_continuous(mu0, gp, ep)
    % Continuous optimum n_p (may be non-integer):
    % n_p = (1 - g_p) / ( sqrt(mu_0 g_p / (1 - e_p)) - g_p )
    np = (1 - gp) ./ (sqrt(mu0 .* gp./(1 - ep)) - gp);
end

function mu0_min = mu0_min_optimum(gp, ep)
    % Lower bound on mu0 for optimum to have positive biomass:
    % mu_0 > ( sqrt(e_p (1 - g_p)) + sqrt(g_p (1 - e_p)) )^2
    if gp>1 || ep>1
        error('The function requires the input g_p<1, e_p<1.');
    end
    mu0_min = (sqrt(ep*(1 - gp)) + sqrt(gp*(1 - ep))).^2;
end

function mu0_min = mu0_min_SS(NT, gp, ep)
    % Lower bound on mu0 for SS to have positive biomass
    % mu_0 > (1 + (N_T-1) g_p)(1 + (N_T-1) e_p) / N_T
    mu0_min = (1 + (NT - 1)*gp) .* (1 + (NT - 1)*ep) ./ NT;
end


function CustomColormap = colormapGeneration(MIN_dM, MAX_dM)
colorstep = 1e-4;
ncolor = floor(MIN_dM/colorstep):ceil(MAX_dM/colorstep);
izero = find(ncolor == 0);
nn = izero - 1;
np = length(ncolor) - izero;
CustomColormap = NaN*ones(length(ncolor),3);
colorn = [0.20 0.20 0.90];
colorp = [1.00 0.90 0.10];

alpha_color = 0.00;

for i = 1:length(ncolor)
    switch sign(i - izero)
        case -1
            CustomColormap(i,:) = colorn + (1-alpha_color) * (1 - colorn) * i/nn;
        case 0
            CustomColormap(i,:) = [1 1 1];
        case 1
            CustomColormap(i,:) = (1-alpha_color) + alpha_color * colorp - (1-alpha_color) * (1 - colorp) * (i-izero)/np;
    end
end
CustomColormap = CustomColormap*0.95;
end

