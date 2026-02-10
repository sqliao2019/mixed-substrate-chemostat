clear;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'data');

dataFile = fullfile(dataDir,'data_Fig3.mat');

if ~exist(dataFile,'file')
    error(['Required data file not found: %s\n' ...
           'Please run Fig3_1_run.m first, or download the precomputed data.'], ...
           dataFile);
end
%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%% 
load(dataFile)
scan = data_Fig3.scan;

%% Plot panel b
figure; hold on;
plot(scan.eta_values, data_Fig3.b.M(1,:),'--k');
plot(scan.eta_values, data_Fig3.b.M(2,:),'--k');
plot(scan.eta_values, data_Fig3.b.M(1,:)+data_Fig3.b.M(2,:),'--r');
plot(scan.eta_values, data_Fig3.b.M(3,:),'b');
title('Fig3b');
xlabel('Fraction of input substrate 1');
ylabel('Biomass');
axis([0 1 0 8])
box on;

%% Compute the shared colorbar for panel c and d
%%% M1ave and M2ave have slight difference due to asymmetry of the first
%%% condition. To avoid discussing this feature in the main text, Fig3c
%%% plots 2*M1ave instead of M1ave + M2ave
Md_ave = 2 * data_Fig3.c.Md1ave;
Ms_ave = data_Fig3.d.Msave;
winner = sign(Md_ave - Ms_ave);
winner(abs(Md_ave - Ms_ave)<0.01) = 0;

val = [Md_ave(:);Ms_ave(:)];
cmin = min(val);
cmax = max(val);

%% Plot panel c
figure; hold on;
pcolor(scan.L_values, log10(scan.F_values), 2 * data_Fig3.c.Md1ave); shading flat; % shading flat to display the discontinuous drop 
plot(0.25, log10(2), 'o', 'MarkerEdgeColor','None', 'MarkerFaceColor','r');
plot(0.20, log10(1e-3), 'o', 'MarkerEdgeColor','None', 'MarkerFaceColor','r');
plot(0.40, log10(1e-3), 'o', 'MarkerEdgeColor','None', 'MarkerFaceColor','r');
colorbar;
clim([cmin cmax]);
title('Fig3c');
xlabel('Amplitude L');
ylabel('Frequency F');
box on;

%% Plot panel d
figure; hold on;
pcolor(scan.L_values, log10(scan.F_values), data_Fig3.d.Msave); shading flat;
plot(0.25, log10(2), 'o', 'MarkerEdgeColor','None', 'MarkerFaceColor','r');
plot(0.20, log10(1e-3), 'o', 'MarkerEdgeColor','None', 'MarkerFaceColor','r');
plot(0.40, log10(1e-3), 'o', 'MarkerEdgeColor','None', 'MarkerFaceColor','r');
colorbar;
clim([cmin cmax]);
title('Fig3d');
xlabel('Amplitude L');
ylabel('Frequency F');
box on;

%% Plot panel e
figure; hold on;

pcolor(scan.L_values, log10(scan.F_values), winner); shading flat; 
plot(0.25, log10(2), 'o', 'MarkerEdgeColor','None', 'MarkerFaceColor','r');
plot(0.20, log10(1e-3), 'o', 'MarkerEdgeColor','None', 'MarkerFaceColor','r');
plot(0.40, log10(1e-3), 'o', 'MarkerEdgeColor','None', 'MarkerFaceColor','r');
title('Fig3e')
xlabel('Amplitude L');
ylabel('Frequency F');
clim([-1 1])
colorp = [0.20 0.20 0.90];
colorn = [1.00 0.90 0.00];
colorz = [1.00 1.00 1.00];

colormap([colorp;colorz;colorn]*0.95);

colorbar;
box on;

%% Plot panel f
figure; hold on;
plot(data_Fig3.f.T,data_Fig3.f.Y(:,7),'b-');
plot(data_Fig3.f.T,data_Fig3.f.Y(:,6),'k-');
plot(data_Fig3.f.T,data_Fig3.f.Y(:,5)+data_Fig3.f.Y(:,6),'r-');  
title('Fig3f');
xlabel('Time');
ylabel('Biomass');
ylim([0 8])
xlim([0 10])
box on;

%% Plot panel g
figure; hold on;

a = 20; b = 500; r = 1/120;
t  = data_Fig3.g.T;  x = timeWarp(t, b, a, r);

plot(x,data_Fig3.g.Y(:,7),'b-');
plot(x,data_Fig3.g.Y(:,6),'k-');
plot(x,data_Fig3.g.Y(:,5)+data_Fig3.g.Y(:,6),'r-');  
title('Fig3g');
xlabel('Time');
ylabel('Biomass');
marks = [0 5 10 15 20 500 505 510 515 520];
xticks(timeWarp(marks, b, a, r)); xticklabels(string(marks));
xlim([timeWarp(marks(1), b, a, r) timeWarp(marks(end), b, a, r)])
ylim([0 8]); 
box on;

%% Plot panel h
figure; hold on;

a = 20; b = 500; r = 1/120;
t  = data_Fig3.h.T;  x = timeWarp(t, b, a, r);

plot(x,data_Fig3.h.Y(:,7),'b-');
plot(x,data_Fig3.h.Y(:,6),'k-');
plot(x,data_Fig3.h.Y(:,5)+data_Fig3.h.Y(:,6),'r-');  
title('Fig3h');
xlabel('Time');
ylabel('Biomass');
marks = [0 5 10 15 20 500 505 510 515 520];
xticks(timeWarp(marks, b, a, r)); xticklabels(string(marks));
xlim([timeWarp(marks(1), b, a, r) timeWarp(marks(end), b, a, r)])
ylim([0 8]); 
box on;

%%
% time-warp helper (piecewise linear compress on [a,b] by factor r)
function x = timeWarp(t, half, a, r)
%TIMEWARPPERIODIC  Periodic time warp per half-period.
%  - Keep first 'a' units in each half uncompressed
%  - Compress the remaining (half - a) by factor r
%  - Works for any t vector, no loops

    a = min(a, half);                         % clamp
    n = floor(t ./ half);                     % which half are we in?
    phase = t - n .* half;                    % time within current half
    
    compressed_per_half = a + r * (half - a);
    base_x = n .* compressed_per_half;
    x = base_x + min(phase, a) + r * max(phase - a, 0);
end