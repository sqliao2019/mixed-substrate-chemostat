clear; clc;

baseDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(baseDir,'data');
if ~exist(dataDir,'dir'); mkdir(dataDir); end

%%%%%%%%%%%%%%%%%%%%%%% Setting %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%% Figure-specific setting (Fig3)

% Substrate Uptake Function  
alpha_value = 4.0; Km_value = 2.0;
p0.f1 = @(x) alpha_value .* x./ (x + Km_value);
p0.f2 = @(x) alpha_value .* x./ (x + Km_value);

% Input
St = 10; % total input
Sc = 5;  % averaged input substrate 1

%%%%%% ODE option
opt = odeset('NonNegative',1:7,'RelTol',1e-6,'AbsTol',1e-9);

%%%%%% Other Parameters
%%% Parameters are involved in struct p.

p0.dilution = 2.0; % Dilution rate gamma in main text/SI
p0.Ad1 = 1.0; p0.md1 = 0.0; % DOL 1
p0.Ad2 = 1.0; p0.md2 = 0.0; % DOL 2

% SS
p0.As = 1.0; p0.ms = 0.0; 
p0.sig1 = 0.5; 
p0.sig2 = 0.5;

p0.cutoff = 1e-5; % below which the strain is completely washed out.

%%%%%% helper
init = @(S1i, S2i, M1i, M2i, Msi) [S1i, S2i, S1i, S2i, M1i, M2i, Msi]; % initial condition 

%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%% 
%% Compute the results for panel b
scan.eta_values = 0:0.01:1;
data_Fig3.b.M = NaN* ones(3, numel(scan.eta_values));
S10v = St * scan.eta_values;
S20v = St * (1-scan.eta_values);

for i = 1:numel(scan.eta_values)
    [d1, d2, ss] = solveGAS(p0, S10v(i), S20v(i));
    data_Fig3.b.M(:,i) = [d1(end), d2(end), ss(end)];
end

%% Compute the results for panel c-d
scan.F_values = logspace(1, -3, 201);
scan.L_values = 0:0.002:min(Sc/St, 1-Sc/St);

n1 = numel(scan.F_values);
n2 = numel(scan.L_values);
data_Fig3.c.Md1ave = NaN * ones(n1,n2);
data_Fig3.c.Md2ave = NaN * ones(n1,n2);
data_Fig3.d.Msave = NaN * ones(n1,n2);

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
        
        data_Fig3.c.Md1ave(i,j) = 1/2 * (Md1(1:end-1) + Md1(2:end))' * diff(T) / Tperiod; % time average of M1
        data_Fig3.c.Md2ave(i,j) = 1/2 * (Md2(1:end-1) + Md2(2:end))' * diff(T) / Tperiod; % time average of M2
        data_Fig3.d.Msave(i,j) = 1/2 * (Ms(1:end-1) + Ms(2:end))' * diff(T) / Tperiod; % time average of Ms
    end
end

%% Compute the result for panel f
F = 2; L = 0.25;

Tperiod = 1/F;

S10a = Sc + St * L; 
S20a = St - Sc - St * L; 
S10b = Sc - St * L; 
S20b = St - Sc + St * L;

y0 = init(Sc, St - Sc, 0.1, 0.1, 0.2);

Tf = 0;
Yf = y0;
t0 = 0;

Tspan_half = makeTspanHalf(Tperiod/2, 0.001, 12, 20);

for pn = 1:20
    [t1,y1] = ode15s(@(t,y)model_general(t,y,p0,S10a,S20a), t0 + Tspan_half, y0, opt);
    y0 = y1(end,:);
    [t2,y2] = ode15s(@(t,y)model_general(t,y,p0,S10b,S20b), t0 + Tspan_half + Tperiod/2, y0, opt);
    y0 = y2(end,:);
    t0 = t2(end);
    Tf = [Tf(1:end-1); t1(1:end-1); t2];
    Yf = [Yf(1:end-1,:); y1(1:end-1,:); y2];  
end
Yf(:,5:7) = Yf(:,5:7) .* (Yf(:,5:7) > p0.cutoff);
data_Fig3.f.T = Tf;
data_Fig3.f.Y = Yf;

%% Compute the result for panel g
F = 1e-3; L = 0.2;

Tperiod = 1/F;

S10a = Sc + St * L; 
S20a = St - Sc - St * L; 
S10b = Sc - St * L; 
S20b = St - Sc + St * L;

y0 = init(Sc, St - Sc, 0.1, 0.1, 0.2);

Tg = 0;
Yg = y0;
t0 = 0;

Tspan_half = makeTspanHalf(Tperiod/2, 0.1, 12, 20);

for pn = 1:10
    [t1,y1] = ode15s(@(t,y)model_general(t,y,p0,S10a,S20a), t0 + Tspan_half, y0, opt);
    y0 = y1(end,:);
    [t2,y2] = ode15s(@(t,y)model_general(t,y,p0,S10b,S20b), t0 + Tspan_half + Tperiod/2, y0, opt);
    y0 = y2(end,:);
    t0 = t2(end);
    Tg = [Tg(1:end-1); t1(1:end-1); t2];
    Yg = [Yg(1:end-1,:); y1(1:end-1,:); y2];  
end
Yg(:,5:7) = Yg(:,5:7) .* (Yg(:,5:7) > p0.cutoff);
data_Fig3.g.T = Tg;
data_Fig3.g.Y = Yg;

%% Compute the result for panel h
F = 1e-3; L = 0.4;

Tperiod = 1/F;

S10a = Sc + St * L; 
S20a = St - Sc - St * L; 
S10b = Sc - St * L; 
S20b = St - Sc + St * L;

y0 = init(Sc, St - Sc, 0.1, 0.1, 0.2);

Th = 0;
Yh = y0;
t0 = 0;

Tspan_half = makeTspanHalf(Tperiod/2, 0.1, 12, 20);

for pn = 1:10
    [t1,y1] = ode15s(@(t,y)model_general(t,y,p0,S10a,S20a), t0 + Tspan_half, y0, opt);
    y0 = y1(end,:);
    [t2,y2] = ode15s(@(t,y)model_general(t,y,p0,S10b,S20b), t0 + Tspan_half + Tperiod/2, y0, opt);
    y0 = y2(end,:);
    t0 = t2(end);
    Th = [Th(1:end-1); t1(1:end-1); t2];
    Yh = [Yh(1:end-1,:); y1(1:end-1,:); y2];  
end
Yh(:,5:7) = Yh(:,5:7) .* (Yh(:,5:7) > p0.cutoff);
data_Fig3.h.T = Th;
data_Fig3.h.Y = Yh;



%% Save data
data_Fig3.scan = scan;
save(fullfile(dataDir,'data_Fig3.mat'),'data_Fig3');

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

function Tspan = makeTspanHalf(halfPeriod, fineStep, coarseStep, maxFineRange)
    % generate hybrid time span with dense sampling near start
    fineEnd = min(maxFineRange, halfPeriod);
    Tspan = [0:fineStep:fineEnd-fineStep, ...
                 fineEnd:coarseStep:halfPeriod-coarseStep, ...
                 halfPeriod];
end