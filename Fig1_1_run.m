clear; clc;

baseDir = fileparts(mfilename('fullpath'));
autoDir = fullfile(baseDir,'autogenFig1');
dataDir = fullfile(baseDir,'data');
if ~exist(autoDir,'dir'); mkdir(autoDir); end
if ~exist(dataDir,'dir'); mkdir(dataDir); end

%%%%%%%%%%%%%%%%%%%%%%% Setting %%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%% Parameters
%%% Parameters are involved in struct p.

p0.dilution = 2.5; % Dilution rate gamma in main text/SI
p0.Ad1 = 1.0; p0.md1 = 0.5; % DOL 1
p0.Ad2 = 1.0; p0.md2 = 1.0; % DOL 2

% SS
p0.As = 1.5;   p0.ms = 1.0; 
p0.sig1 = 1.0; p0.sig2 = 1.0;

% Substrate Uptake Function
alpha1 = 4; K1 = 2;
alpha2 = 6; K2 = 3;

p0.f1 = @(x) alpha1 .*x ./(x+K1);
p0.f2 = @(x) alpha2 .*x ./(x+K2);

%%%%% Initial Substrate Input
S10_values = 0:0.01:10;
S20_values = 0:0.01:10;

%%%%% Scanned factors, S1 and S2
S1_values = S10_values;
S2_values = S20_values;

%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%% 
%% Compute results for Figure 1e
data_Fig1.e.results_Lamd1 = (p0.f1(S1_values) - p0.md1)/p0.Ad1;
data_Fig1.e.results_Lamd2 = (p0.f2(S2_values) - p0.md2)/p0.Ad2;

syms x1 x2
eqn1 = [(p0.f1(x1) - p0.md1)/p0.Ad1 == p0.dilution, ...
        (p0.f2(x2) - p0.md2)/p0.Ad2 == p0.dilution];

sol_1 = solve(eqn1, [x1, x2]);
data_Fig1.e.Cd1 = double(sol_1.x1);
data_Fig1.e.Cd2 = double(sol_1.x2);

%% Compute results for Figure 1f
data_Fig1.f.results_Lams = (p0.sig1 * p0.f1(S1_values) + p0.sig2 * (p0.f2(S2_values))' - p0.ms) / p0.As; % f1 + f2'

%% Compute results for Figure 1g-h
% obtain the expression of solution in symbols
syms y [7 1] 
syms varS10 varS20
sol = solve(model_general(0, y, p0, varS10, varS20)==0, y);

% solution for DOL1
steadyFun.DOL1 = matlabFunction(...
    unique([sol.y1, sol.y5], 'rows', 'stable'), ...
    'Vars', {varS10}, ...    
    'File',  fullfile(autoDir,'DOL1_ss'));
% solution for DOL2
steadyFun.DOL2 = matlabFunction(...
    unique([sol.y2, sol.y6], 'rows', 'stable'), ...
    'Vars', {varS20}, ...    
    'File',  fullfile(autoDir,'DOL2_ss'));
% solution for SS
steadyFun.SS   = matlabFunction(...
    unique([sol.y3, sol.y4, sol.y7], 'rows', 'stable'), ...
    'Vars', {varS10, varS20}, ...    
    'File',  fullfile(autoDir,'SS_ss'));

addpath(autoDir);
%%
n1 = numel(S10_values);
n2 = numel(S20_values);
data_Fig1.g.results_Md1 = NaN(1, n1);
data_Fig1.g.results_Md2 = NaN(1, n2);
data_Fig1.h.results_Ms = NaN(n2, n1);

for i = 1: n1
    S10 = S10_values(i);
    sol_row = selectSteadyState(steadyFun.DOL1(S10));
    data_Fig1.g.results_Md1(i) = sol_row(end);
end

for i = 1: n2
    S20 = S20_values(i);  
    sol_row = selectSteadyState(steadyFun.DOL2(S20));
    data_Fig1.g.results_Md2(i) = sol_row(end);
end

for i = 1: n1
    S10 = S10_values(i);
    for j = 1: n2
        S20 = S20_values(j);
        sol_row = selectSteadyState(steadyFun.SS(S10,S20));
        data_Fig1.h.results_Ms(j, i) = sol_row(end);
    end
end

if isfolder(autoDir)
    try
        % recursive remove (deletes all files/subfolders)
        rmpath(autoDir);
        rmdir(autoDir, 's');
    catch ME
        warning('Could not delete folder %s:\n%s', autoDir, ME.message);
    end
end

%% Save data
data_Fig1.scan.S1_values = S1_values;
data_Fig1.scan.S2_values = S2_values;
data_Fig1.scan.S10_values = S10_values;
data_Fig1.scan.S20_values = S20_values;
data_Fig1.p0 = p0;
save(fullfile(dataDir,'data_Fig1.mat'),'data_Fig1');


%%
function sol_row = selectSteadyState(sol)
    % sol: candidate steady states, each row is one solution, last column = biomass
    mask = all(sol >= 0, 2);
    sol = sol(mask, :);
    if isempty(sol)
        error('No non-negative steady state found (expected at least a trivial one).');
    end
    sol = unique(sol, 'rows', 'stable');
    if size(sol,1) > 2
        warning('Multiple non-trivial steady states found. Using the one with maximum biomass.');
    end
    [~, idx] = max(sol(:,end));
    sol_row = sol(idx,:);
end