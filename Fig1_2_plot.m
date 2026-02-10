clear;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'data');

dataFile = fullfile(dataDir,'data_Fig1.mat');

if ~exist(dataFile,'file')
    error(['Required data file not found: %s\n' ...
           'Please run Fig1_1_run.m first, or download the precomputed data.'], ...
           dataFile);
end
%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%% 
load(dataFile,'data_Fig1');
S1_values = data_Fig1.scan.S1_values;
S2_values = data_Fig1.scan.S2_values;
S10_values = data_Fig1.scan.S10_values;
S20_values = data_Fig1.scan.S20_values;
p0 = data_Fig1.p0;

%% Plot Figure 1e
figure; hold on;
plot(S1_values, data_Fig1.e.results_Lamd1);
plot(S2_values, data_Fig1.e.results_Lamd2);
plot([0 10],[p0.dilution p0.dilution]);
plot(data_Fig1.e.Cd1, p0.dilution,'o');
plot(data_Fig1.e.Cd2, p0.dilution,'o');
plot([data_Fig1.e.Cd1, data_Fig1.e.Cd1], [0, p0.dilution]);
plot([data_Fig1.e.Cd2, data_Fig1.e.Cd2], [0, p0.dilution]);
ylim([0,4]);
box on;
xticks(0)
yticks(0)
title('Fig. 1e')
xlabel('Substrate concentration');
ylabel('Cell growth rate');

%% Plot Figure 1f
figure; hold on;
contourf(S1_values, S2_values, data_Fig1.f.results_Lams, -1:0.5:5, 'LineStyle','none'); 
contour(S1_values, S2_values, data_Fig1.f.results_Lams, [p0.dilution, p0.dilution], 'k');
cmap1 = linspace(1.0,0.2,12)'*[1 1 0]+[0 0 1];
colormap(cmap1)
clim([-1 5])
colorbar;
box on;
title('Fig. 1f')

%% Plot Figure 1g
figure; hold on;
plot(S10_values, data_Fig1.g.results_Md1);
plot(S20_values, data_Fig1.g.results_Md2);
box on;
xticks(0)
yticks(0)
title('Fig. 1g')

%% Plot Figure 1h
figure; hold on;
Ms = data_Fig1.h.results_Ms;
Ms(Ms == 0) = NaN;
contourf(S10_values, S20_values, Ms, 0:0.5:8.5, 'LineStyle','none');
cmap2 = linspace(1.0,0.0,18)'*[1 1 0]+[0 0 1];
cmap2(1,:) = [];
colormap(cmap2)
clim([0 8.5])
colorbar;
box on;
contour(S10_values, S20_values, data_Fig1.f.results_Lams, [p0.dilution, p0.dilution], 'k');
title('Fig. 1h')
xlim([min(S10_values) max(S10_values)])
ylim([min(S20_values) max(S20_values)])


