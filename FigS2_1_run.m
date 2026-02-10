clear; clc;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'data');
if ~exist(dataDir,'dir'); mkdir(dataDir); end
%%%%%%%%%%%%%%%%%%%%%%% Setting %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%% Figure-specific setting (FigS2)

% Substrate Uptake Function  
alpha_value = 4.0; Km_value = 2.0;
p0.f1 = @(x) alpha_value .* x./ (x + Km_value);
p0.f2 = @(x) alpha_value .* x./ (x + Km_value);

% Input
St = 10; % total input
Sc = 3;  % averaged input substrate 1

%%%%%% ODE option
opt = odeset('NonNegative',1:7,'RelTol',1e-6,'AbsTol',1e-9);

%%%%%% Other Parameters
%%% Parameters are involved in struct p.

p0.dilution = 2.0; % Dilution rate gamma in main text/SI
p0.Ad1 = 1.0; p0.md1 = 0.0; % DOL 1
p0.Ad2 = 1.0; p0.md2 = 0.0; % DOL 2

% SS
p0.As = 1.0; p0.ms = 0.0; 
sig_crit = computeCritSigma(p0, Sc, St-Sc);
p0.sig1 = sig_crit; 
p0.sig2 = sig_crit;

p0.cutoff = 1e-5; % below which the strain is completely washed out.

%%%%%% helper
init = @(S1i, S2i, M1i, M2i, Msi) [S1i, S2i, S1i, S2i, M1i, M2i, Msi]; % initial condition 

%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%% 
%% Compute the results for panel b
scan.eta_values = 0:0.01:1;
data_FigS2.b.M = NaN* ones(3, numel(scan.eta_values));
S10v = St * scan.eta_values;
S20v = St * (1-scan.eta_values);

for i = 1:numel(scan.eta_values)
    [d1, d2, ss] = solveGAS(p0, S10v(i), S20v(i));
    data_FigS2.b.M(:,i) = [d1(end), d2(end), ss(end)];
end

%% Compute the results for panel c-e
scan.F_values = logspace(1, -3, 201);
scan.L_values = 0:0.002:min(Sc/St, 1-Sc/St);

n1 = numel(scan.F_values);
n2 = numel(scan.L_values);
data_FigS2.c.Md1ave = NaN * ones(n1,n2);
data_FigS2.c.Md2ave = NaN * ones(n1,n2);
data_FigS2.d.Msave = NaN * ones(n1,n2);

for i = 1:n1
    for j = 1:n2
        %%% Period
        Tperiod = 1.0/scan.F_values(i);
        %%% Substrate input
        % First half of period
        S10a = Sc + St * scan.L_values(j); 
        S20a = St - Sc - St * scan.L_values(j); 
        % Second half of period
        S10b = Sc - St * scan.L_values(j); 
        S20b = St - Sc + St * scan.L_values(j);
        % Initial condition
        y0 = init(Sc, St - Sc, 0.1, 0.1, 0.2);

        nWarm = floor(100/Tperiod + 100);
        for pn = 1:nWarm
            [ta,ya] = ode15s(@(t,y)model_general(t,y,p0,S10a,S20a),[0 Tperiod/2], y0, opt);
            y0 = ya(end,:);
            [tb,yb] = ode15s(@(t,y)model_general(t,y,p0,S10b,S20b),[Tperiod/2 Tperiod], y0, opt);
            y0 = yb(end,:);
        end
        T = [ta(1:end-1);tb];
        assert(all(diff(T) > 0), 'Non-monotone time grid detected.');
        
        Y = [ya(1:end-1,:);yb];
        Y(:,5:7) = Y(:,5:7) .* (Y(:,5:7)>p0.cutoff);
        Md1 = Y(:,5); 
        Md2 = Y(:,6); 
        Ms = Y(:,7);
        
        data_FigS2.c.Md1ave(i,j) = 1/2 * (Md1(1:end-1) + Md1(2:end))' * diff(T) / Tperiod; % time average of M1
        data_FigS2.c.Md2ave(i,j) = 1/2 * (Md2(1:end-1) + Md2(2:end))' * diff(T) / Tperiod; % time average of M2
        data_FigS2.d.Msave(i,j) = 1/2 * (Ms(1:end-1) + Ms(2:end))' * diff(T) / Tperiod; % time average of Ms
    end
end


%% Save data
data_FigS2.scan = scan;
save(fullfile(dataDir,'data_FigS2.mat'),'data_FigS2');

%%
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

function sig_crit = computeCritSigma(p0, S10, S20)
    % Calculate the critical value of sigma that allows Ms = Md at g = 1, 
    % S10 = Sc, S20 = St-Sc.
    syms y [7 1] real
    syms varsig real
    p0.sig1 = varsig;
    p0.sig2 = varsig;
    idx_dol1 = [1 5];
    idx_dol2 = [2 6];
    idx_ss = [3 4 7];
    eqns = model_general(0, y, p0, S10, S20)==0;
    eqns_dol1 = eqns(idx_dol1);
    eqns_dol2 = eqns(idx_dol2);
    eqns_ss = eqns(idx_ss);
    
    sol_dol1 = solve(eqns_dol1, y(idx_dol1));
    sol_dol1_row = selectSteadyState(double([sol_dol1.y1, sol_dol1.y5]));
    sol_dol2 = solve(eqns_dol2, y(idx_dol2));
    sol_dol2_row = selectSteadyState(double([sol_dol2.y2, sol_dol2.y6]));

    Md = sol_dol1_row(end) + sol_dol2_row(end);

    if Md <= 0
        error('Non-positive Md. For comparison, please select parameters allowing positive Md.');
    end

    sol_ss = solve([eqns_ss; y(7) == Md; varsig>=0; varsig<=1; y(idx_ss)>=0], [y(idx_ss); varsig]);
    sig_crit = double(sol_ss.varsig);
end