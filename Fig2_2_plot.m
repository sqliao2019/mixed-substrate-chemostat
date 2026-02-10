clear;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'data');

dataFile = fullfile(dataDir,'data_Fig2.mat');

if ~exist(dataFile,'file')
    error(['Required data file not found: %s\n' ...
           'Please run Fig2_1_run.m first, or download the precomputed data.'], ...
           dataFile);
end
%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%% 
load(dataFile)
scan = data_Fig2.scan;

%% Plot Figure 2a
figure; hold on;
MIN_dM = min(data_Fig2.a.DeltaM, [], 'all');
MAX_dM = max(data_Fig2.a.DeltaM, [], 'all');
CustomColormap = colormapGeneration(MIN_dM, MAX_dM);
pcolor(scan.g_values, scan.sig_values, data_Fig2.a.DeltaM); shading interp;
clim([MIN_dM, MAX_dM]);
colormap(CustomColormap);
plot(data_Fig2.a.gEBB, scan.sig_values, 'k');
plot(data_Fig2.a.gRGB, scan.sig_values, 'k-');
title('Fig. 2a');
xlabel('Relative metabolic cost g');
ylabel('Relative uptake efficiency \sigma');
axis([1 2 0 1]);
box on;

%% Plot Figure 2b
figure; hold on;
MIN_dM = min(data_Fig2.a.DeltaM, [], 'all'); %% Using the same colormap as Fig.2a for comparison
MAX_dM = max(data_Fig2.a.DeltaM, [], 'all');
CustomColormap = colormapGeneration(MIN_dM, MAX_dM);
pcolor(scan.g_values, scan.sig_values, data_Fig2.b.DeltaM); shading interp;
clim([MIN_dM, MAX_dM]);
colormap(CustomColormap);
% plot(data_Fig2.b.gEBB, scan.sig_values, 'k'); % meaningless due to negative biomass.
plot(data_Fig2.b.gSWB, scan.sig_values, 'k');
plot(data_Fig2.a.gRGB, scan.sig_values, 'k-');
title('Fig. 2b');
xlabel('Relative metabolic cost g');
ylabel('Relative uptake efficiency \sigma');
axis([1 2 0 1]);
box on;

%% Plot Figure 2c
figure; hold on;
MIN_dM = min(data_Fig2.c.DeltaM, [], 'all');
MAX_dM = max(data_Fig2.c.DeltaM, [], 'all');
CustomColormap = colormapGeneration(MIN_dM, MAX_dM);
pcolor(scan.mu_values, scan.h_values, data_Fig2.c.DeltaM); shading interp; % using g as y axis (note that here g = 1+h, sigma = 1-h) 
clim([MIN_dM, MAX_dM]);
colormap(CustomColormap);
plot(data_Fig2.c.muDWB, scan.h_values, 'k');
plot(data_Fig2.c.muSWB, scan.h_values, 'k');
plot(data_Fig2.c.muEBB, scan.h_values, 'k');
plot(scan.mu_values, data_Fig2.c.hRGB, 'k');
title('Fig. 2c');
xlabel('Normalized total substrate \mu');
ylabel('\sigma/g');
axis([1 10 0 1])
box on;

%% Plot Figure 2d
figure; hold on;
plot(scan.eta_values, data_Fig2.d.Md(1,:),'r');
plot(scan.eta_values, data_Fig2.d.Ms,'b');
title('Fig. 2d');
xlabel('Fraction of input substrate 1');
ylabel('Biomass');
axis([0 1 0 5])
box on;

%% Plot Figure 2e
figure; hold on;
plot(scan.eta_values, data_Fig2.e.Md(1,:),'r');
plot(scan.eta_values, data_Fig2.e.Ms,'b');
title('Fig. 2e');
xlabel('Fraction of input substrate 1');
ylabel('Biomass');
axis([0 1 0 5])
box on;

%% Plot Figure 2f
figure; hold on;
plot(scan.eta_values, data_Fig2.f.Md(1,:),'r');
plot(scan.eta_values, data_Fig2.f.Ms,'b');
title('Fig. 2f');
xlabel('Fraction of input substrate 1');
ylabel('Biomass');
axis([0 1 0 5])
box on;

%%
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
