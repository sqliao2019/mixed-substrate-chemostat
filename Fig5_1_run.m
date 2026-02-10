clear; clc;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'data');
if ~exist(dataDir,'dir'); mkdir(dataDir); end

%%%%%%%%%%%%%%%%%%%%%%% Setting %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%% Shared Parameters
%%% Parameters are involved in struct p.

p0.dilution = 1.0; % Dilution rate
p0.Ad1 = 1.0;    p0.md1 = 0.0;    % DOL 1
p0.Ad2 = p0.Ad1; p0.md2 = p0.md1; % DOL 2
p0.Ad3 = p0.Ad1; p0.md3 = p0.md1; % DOL 3

% SS
p0.ms = 0.0;
% p.As, p.sig1, p.sig2 and p.sig3 are initialzed at the block of each panel

p0.cutoff = 1e-5; % below which the strain is completely washed out.

%%%%%% Panel-specific parameters
% Substrate Uptake Function
p0.alpha1 = 2.5; p0.K1 = 1.0;
p0.alpha2 = 2.5; p0.K2 = 1.0;
p0.alpha3 = 2.5; p0.K3 = 1.0;
p0.f1 = @(x) p0.alpha1.*x ./ (x+p0.K1);
p0.f2 = @(x) p0.alpha2.*x ./ (x+p0.K2);
p0.f3 = @(x) p0.alpha3.*x ./ (x+p0.K3);

%%%%% input
St = 3;

%%%%%% Substrate composition (mapping R1+R2+R3=1 to a 2D equilateral triangle)
m = 1001;
n = 501;
% m and n should be both odd;
m = m + 1 - mod(m,2);
n = n + 1 - mod(n,2);

Xv = linspace(0,1,m);
Yv = linspace(0,sqrt(3)/2,n);
grid.X = Xv;
grid.Y = Yv;

%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%% 


%% Calculate results for Fig5a-d (For Fig5e-i, the results are analytically computed in Fig5_2_plot.m)
gScan = [1.0, 1.2, 2.0];
w = length(gScan);
data_Fig5.Md = NaN * ones(n,m,w);
data_Fig5.Ms = NaN * ones(n,m,w);

z0 = [];
for k = 1:w
    p0.As = p0.Ad1 * gScan(k);
    p0.sig1 = 1;
    p0.sig2 = 1;
    p0.sig3 = 1;
    for j = 1:m
        for i = 1:n
            % Coordinates on substrate composition triangle
            
            Xc = Xv(j); Yc = Yv(i);
            if Yc > sqrt(3) * min(Xc, 1-Xc)
                continue;
            end
            %[i j]
            d1 = Yc;
            d2 = abs(sqrt(3)*Xc/2-Yc/2); 
            R1 = d1 * 2/sqrt(3); % rescale R1_max to 1
            R2 = d2 * 2/sqrt(3);
            R3 = 1 - R1 - R2;
            S10 = R1 * St;
            S20 = R2 * St;
            S30 = R3 * St;
            %[R1 R2 R3]
            sol = solve_3dol_ss_threeSubstrate(p0, S10, S20, S30);

            data_Fig5.Md(i,j,k) = sum(sol(7:9));
            data_Fig5.Ms(i,j,k) = sol(10);
        end
    end
end

%%
data_Fig5.grid = grid;
save(fullfile(dataDir,'data_Fig5.mat'), 'data_Fig5');


