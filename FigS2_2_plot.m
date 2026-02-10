clear;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'data');

dataFile = fullfile(dataDir,'data_FigS2.mat');

if ~exist(dataFile,'file')
    error(['Required data file not found: %s\n' ...
           'Please run FigS2_1_run.m first, or download the precomputed data.'], ...
           dataFile);
end
%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%% 
load(dataFile)
scan = data_FigS2.scan;

%% Plot panel b
figure; hold on;
plot(scan.eta_values, data_FigS2.b.M(1,:),'--k');
plot(scan.eta_values, data_FigS2.b.M(2,:),'--k');
plot(scan.eta_values, data_FigS2.b.M(1,:)+data_FigS2.b.M(2,:),'--r');
plot(scan.eta_values, data_FigS2.b.M(3,:),'b');
title('FigS2b');
xlabel('Fraction of input substrate 1');
ylabel('Biomass');
axis([0 1 0 8])
box on;

%% Compute the shared colorbar for panel c and d
Md_ave = data_FigS2.c.Md1ave + data_FigS2.c.Md2ave;
Ms_ave = data_FigS2.d.Msave;
winner = sign(Md_ave - Ms_ave);
winner(abs(Md_ave - Ms_ave)<0.01) = 0;

val = [Md_ave(:);Ms_ave(:)];
cmin = min(val);
cmax = max(val);

%% Plot panel c
figure; hold on;
pcolor(scan.L_values, log10(scan.F_values), Md_ave); shading flat; % shading flat to display the discontinuous drop 
colorbar;
clim([cmin cmax]);
box on;

%% Plot panel d
figure; hold on;
pcolor(scan.L_values, log10(scan.F_values), Ms_ave); shading flat;
colorbar;
clim([cmin cmax]);
box on;

%% Plot panel e
figure; hold on;
pcolor(scan.L_values, log10(scan.F_values), winner); shading flat; 
clim([-1 1])
colorp = [0.20 0.20 0.90];
colorn = [1.00 0.90 0.00];
colorz = [1.00 1.00 1.00];

colormap([colorp;colorz;colorn]*0.95);

colorbar;
box on;