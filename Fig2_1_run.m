clear; clc;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'data');
if ~exist(dataDir,'dir'); mkdir(dataDir); end

%%%%%%%%%%%%%%%%%%%%%%% Setting %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%% Shared Parameters
%%% Parameters are involved in struct p.

p0.dilution = 2.0; % Dilution rate gamma in main text/SI
p0.Ad1 = 1.0;    p0.md1 = 0.0;   % DOL 1
p0.Ad2 = p0.Ad1; p0.md2 = p0.md1; % DOL 2

% SS
p0.ms = 0.0;
% p.As, p.sig1, p.sig2 are initialized at the block of each panel

%%%%%% Panel-specific parameters

%%% Model with linear uptake rate (used in panel a-d)
p1 = p0;
% Substrate Uptake Function
p1.k_value = 2.0;
% closed form of uptake rate (not directly used here):
% p1.f1 = @(x) p1.k_value.*x;
% p1.f2 = @(x) p1.k_value.*x;

%%% Model with Monod uptake rate (used in panel e)
p2 = p0;
alpha_value = 3; Km_value = 0.5;
p2.f1 = @(x) alpha_value.*x ./ (x+Km_value);
p2.f2 = @(x) alpha_value.*x ./ (x+Km_value);

%%% Model with quadratic uptake rate (used in panel f)
p3 = p0;
q_value = 2.0;
p3.f1 = @(x) q_value.*x.^2;
p3.f2 = @(x) q_value.*x.^2;

%%%%%% Scanned factors, g, sigma, mu and eta
scan.g_values = 1:0.01:2.2;   % relative metabolic cost
scan.sig_values = 0:0.01:1;   % relative substrate uptake efficiency
scan.mu_values = 1:0.1:10;    % total normalized input
scan.eta_values = 0:0.01:1;   % fraction of substrate 1 in input
scan.h_values = 0:0.01:1;     % Parametric coordinate along g+sigma = 2
scan.g_cases = [1, 1.2, 2.0]; % three values of g used in 2d-f

% mu is defined consistently as (S10+S20)/C_d, with C_d = 1 in all panels by parameter choice.
Cd = 1;
S10 = @(mu, eta) Cd * mu .* eta;
S20 = @(mu, eta) Cd * mu .* (1-eta);

%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%% 
%% analytical solution for core-symmetric linear-uptake case
%%% Note: Analytically, dilution/k_value = C_d/U_d
Md = @(mu, eta) p1.dilution ./ p1.k_value .* (max(mu.*eta - 1, 0) + max(mu.*(1-eta) - 1, 0));
Ms = @(g, sig, mu, eta) p1.dilution ./ p1.k_value  .* max(mu./g - 1./sig, 0);

%% Compute results for Fig.2a
%%% Scan g and sigma, with mu and eta fixed
mu_2a = 5;                              % total normalized input
eta_2a = 0.5;                             % fraction of substrate 1 in input

n1 = numel(scan.sig_values);
n2 = numel(scan.g_values);
data_Fig2.a.DeltaM = NaN * ones(n1, n2);

for i = 1:n1
    for j = 1:n2
        g_2a = scan.g_values(j);
        sig_2a = scan.sig_values(i);
        data_Fig2.a.DeltaM(i,j) = Md(mu_2a, eta_2a) - Ms(g_2a, sig_2a, mu_2a, eta_2a);
    end
end

%%% Compute the boundary
data_Fig2.a.gEBB = eqnEBB1(scan.sig_values, mu_2a);
data_Fig2.a.gRGB = eqnRGB(scan.sig_values);
%% Compute results for Fig.2b
%%% Scan g and sigma, with mu and eta fixed
mu_2b = 1.5;                              
eta_2b = 0.5;                             

n1 = numel(scan.sig_values);
n2 = numel(scan.g_values);
data_Fig2.b.DeltaM = NaN * ones(n1, n2);

for i = 1:n1
    for j = 1:n2
        g_2b = scan.g_values(j);
        sig_2b = scan.sig_values(i);
        data_Fig2.b.DeltaM(i,j) = Md(mu_2b, eta_2b) - Ms(g_2b, sig_2b, mu_2b, eta_2b);
    end
end
%%% Compute the boundary
data_Fig2.b.gSWB = eqnSWB1(scan.sig_values, mu_2b);

%% Compute results for Fig.2c
%%% Scan mu and h, with g = 2./(1+h), sigma = 2.*h./(1+h) and eta = 0.5
eta_2c = 0.5;

n1 = numel(scan.h_values);
n2 = numel(scan.mu_values);
data_Fig2.c.DeltaM = NaN * ones(n1, n2);

for i = 1:n1
    for j = 1:n2
        g_2c = 2. /(1 + scan.h_values(i));
        sig_2c = 2.* scan.h_values(i)./(1 + scan.h_values(i));
        mu_2c = scan.mu_values(j);
        data_Fig2.c.DeltaM(i,j) = Md(mu_2c, eta_2c) - Ms(g_2c, sig_2c, mu_2c, eta_2c);
    end
end
%%% Compute the boundary
data_Fig2.c.muEBB = eqnEBB2(scan.h_values);
data_Fig2.c.muSWB = eqnSWB2(scan.h_values);
data_Fig2.c.muDWB = 2 * ones(size(scan.h_values));
data_Fig2.c.hRGB = 0.5 * ones(size(scan.mu_values));
%% Compute results for Fig.2d
%%% Scan eta, with fixed mu and sigma, and g = 1.0, 1.2, 2.0

mu_2d = 5;
sig_2d = 1;

n1 = numel(scan.g_cases);
n2 = numel(scan.eta_values);
data_Fig2.d.Md = NaN * ones(1, n2);
data_Fig2.d.Ms = NaN * ones(n1, n2);

for j = 1:n2
    eta_2d = scan.eta_values(j);
    data_Fig2.d.Md(1,j) = Md(mu_2d, eta_2d);
    for i = 1:n1
        g_2d = scan.g_cases(i);
        data_Fig2.d.Ms(i,j) = Ms(g_2d, sig_2d, mu_2d, eta_2d);
    end
end

%% Compute results for Fig.2e
%%% Monod uptake rate (using parameter set p2)
%%% Scan eta, with fixed mu and sigma, and g = 1.0, 1.2, 2.0
sig_2e = 1;
mu_2e = 5;

p2.sig1 = sig_2e;
p2.sig2 = sig_2e;

n1 = numel(scan.g_cases);
n2 = numel(scan.eta_values);
data_Fig2.e.Md = NaN * ones(n1, n2);
data_Fig2.e.Ms = NaN * ones(n1, n2);

for j = 1:n2
    eta_2e = scan.eta_values(j);   
    S10_2e = S10(mu_2e, eta_2e); 
    S20_2e = S20(mu_2e, eta_2e);   
    for i = 1:n1
        p2.As = scan.g_cases(i) * p2.Ad1;
        [sol_dol1, sol_dol2, sol_ss] = solveGAS(p2, S10_2e, S20_2e);
        data_Fig2.e.Md(i,j) = sol_dol1(end) + sol_dol2(end);
        data_Fig2.e.Ms(i,j) = sol_ss(end);
    end
end

%% Compute results for Fig.2f
%%% Quadratic uptake rate (using parameter set p3)
%%% Scan eta, with fixed mu and sigma, and g = 1.0, 1.2, 2.0
sig_2f = 1;
mu_2f = 5;

p3.sig1 = sig_2f;
p3.sig2 = sig_2f;

n1 = numel(scan.g_cases);
n2 = numel(scan.eta_values);
data_Fig2.f.Md = NaN * ones(n1, n2);
data_Fig2.f.Ms = NaN * ones(n1, n2);

for j = 1:n2
    eta_2f = scan.eta_values(j);   
    S10_2f = S10(mu_2f, eta_2f); 
    S20_2f = S20(mu_2f, eta_2f);  
    for i = 1:n1
        p3.As = scan.g_cases(i) * p3.Ad1;
        [sol_dol1, sol_dol2, sol_ss] = solveGAS(p3, S10_2f, S20_2f);
        data_Fig2.f.Md(i,j) = sol_dol1(end) + sol_dol2(end);
        data_Fig2.f.Ms(i,j) = sol_ss(end);
    end
end

%% Save data

data_Fig2.scan = scan;
save(fullfile(dataDir,'data_Fig2.mat'),'data_Fig2');

%%

function g1 = eqnEBB1(sigma, mu) % equal-biomass boundary on mu = const plane
g1 = 1./(1-(2-1./sigma)./mu);
end

function g2 = eqnSWB1(sigma, mu) % critical SS-washout boundary on mu = const plane
g2 = mu .* sigma;
end

function g3 = eqnRGB(sigma)
g3 = 2 .* sigma;
end

function mu = eqnEBB2(h)  % equal-biomass boundary on g+sigma = 2 plane 
g = 2./(1 + h);
sigma = 2.*h./(1+h);
mu = (2-1./sigma)./(1-1./g);
end

function mu = eqnSWB2(h)  % critical SS-washout boundary on g+sigma = 2 plane 
g = 2./(1 + h);
sigma = 2.*h./(1+h);
mu = g./sigma;
end


function [sol_dol1, sol_dol2, sol_ss] = solveGAS(p, S10, S20)
    % solve the global asymptotically stable solution
    syms y [7 1] 
    sol = solve(model_general(0, y, p, S10, S20)==0, y);
    sol_dol1 = selectSteadyState(double([sol.y1, sol.y5]));
    sol_dol2 = selectSteadyState(double([sol.y2, sol.y6]));
    sol_ss = selectSteadyState(double([sol.y3, sol.y4, sol.y7]));
end

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