clear; clc;

%%%%%%%%%%%%%%%%%%%%%%% Setting %%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%% Shared Parameters
%%% Parameters are involved in struct p.
% Fig.4 parameter note:
% SI uses tildes (~A_d1, ~m_d1, ...) to denote effective parameters in the
% presence of product. In code we keep the original field names (Ad1, md1, 
% ...) for continuity
p0.dilution = 4.0; % Dilution rate gamma in main text/SI
p0.Ad1 = 1.0;    p0.md1 = 1.0;    p0.thd1 = 1.25;    % DOL 1
p0.Ad2 = p0.Ad1; p0.md2 = p0.md1; p0.thd2 = p0.thd1; % DOL 2

% SS
p0.ths = p0.thd1;
% p.As, p.ms, p.sig1, p.sig2 are initialzed before solving ode

% Substrate Uptake Function (Fig4 code expects linear uptake form)
p0.k1 = 2.5;
p0.k2 = 2.5;
p0.f1 = @(x) p0.k1 * x;
p0.f2 = @(x) p0.k2 * x;

% Product interaction factor
nHill = 3;
Kp = 2;
p0.Hill = @(P) P.^nHill ./ (P.^nHill + Kp.^nHill);
p0.phi = @(P, beta) exp(beta .* p0.Hill(P));
% Substrate input
S10 = 5; 
S20 = 5;

%%%%%%%%%%%%%%%%%%%%%%%  Plot %%%%%%%%%%%%%%%%%%%%%%%%% 
%% plot Figure 4c
figure; hold on;
beta_case_symm = [-1.0, -0.5, 0, 0.5, 1.0];
P_scan = 0:0.1:5;
for beta = beta_case_symm
    plot(P_scan, p0.phi(P_scan, beta), 'DisplayName', ['\beta = ',num2str(beta)]);
end
box on;
title('Fig. 4c');
xlabel('Product concentration P');
ylabel('Product regulation factor \phi');
legend('Location','northwest');

%% plot Figure 4d
figure; hold on;
beta_case_symm = [-1.0, -0.5, 0, 0.5, 1.0];
Sd1_scan = 0:0.1:5;
for beta = beta_case_symm
    plot(Sd1_scan, Lam_d1_symm(S10,Sd1_scan, p0, beta), 'DisplayName', ['\beta = ',num2str(beta)]);
end
plot(Sd1_scan, p0.dilution*ones(size(Sd1_scan)), 'DisplayName', 'Dilution rate');
axis([0 5 0 20])
box on;
title('Fig. 4d');
xlabel('Substrate concentration');
ylabel('Steady-state cell growth rate');
legend('Location','northwest');

%% Plot Figure 4e
beta_case_symm = [-1.0, -0.5, 0, 0.5, 1.0];

figure; hold on;
for beta = beta_case_symm
    beta1 = beta;
    beta2 = beta;
    [gebb, sigScan] = computeBoundary_productFeedback(beta1, beta2, p0, S10, S20);
    plot(gebb, sigScan, 'DisplayName', ['\beta = ',num2str(beta)]);
end
axis([1 1.6 0.4 1])
box on;
title('Fig. 4e');
xlabel('Relative metabolic cost g');
ylabel('Relative efficiency of uptake pathway \sigma');
legend('Location','southeast');

%% Plot Figure 4f
beta2 = 0;
figure; hold on;
beta_case_f = [0, 0.5, 1.0, 20.0];
for beta1 = beta_case_f
    [gebb, sigScan] = computeBoundary_productFeedback(beta1, beta2, p0, S10, S20);
    plot(gebb, sigScan, 'DisplayName', ['\beta_1 = ',num2str(beta1)]);
end
axis([1 1.3 0.4 1])
box on;
title('Fig. 4f');
xlabel('Relative metabolic cost g');
ylabel('Relative efficiency of uptake pathway \sigma');
legend('Location','southeast');

%% Plot Figure 4g
beta2 = 0;
beta1_crit = compute_beta1_crit_DOL1washout(p0,S10,S20);
figure; hold on;
beta_case_g = [0, -0.5, beta1_crit];
for beta1 = beta_case_g
    [gebb, sigScan] = computeBoundary_productFeedback(beta1, beta2, p0, S10, S20);
    plot(gebb, sigScan, 'DisplayName', ['\beta_1 = ',num2str(beta1)]);
end
axis([1 1.6 0.4 1])
box on;
title('Fig. 4g');
xlabel('Relative metabolic cost g');
ylabel('Relative efficiency of uptake pathway \sigma');
legend('Location','southeast');

%% Plot Figure 4h
beta2 = 0;
beta1_crit = compute_beta1_crit_DOL1washout(p0,S10,S20);
figure; hold on;
beta_case_h = [beta1_crit, -2.0, -3.0];
for beta1 = beta_case_h
    [gebb, sigScan] = computeBoundary_productFeedback(beta1, beta2, p0, S10, S20);
    plot(gebb, sigScan, 'DisplayName', ['\beta_1 = ',num2str(beta1)]);
end
axis([1 1.6 0.4 1])
box on;
title('Fig. 4h');
xlabel('Relative metabolic cost g');
ylabel('Relative efficiency of uptake pathway \sigma');
legend('Location','southeast');

%%
function lam = Lam_d1_symm(S10, Sd1, param, beta)
    % growth rate with steady-state product-substrate balance map (use DOL 1 as example)
    Ud1 = param.Ad1 + param.md1/param.dilution;
    Md1 = (S10 - Sd1)/Ud1;
    P = 2*param.thd1 * Md1;
    lam = (param.phi(P,beta) .* param.f1(Sd1) - param.md1) / param.Ad1;
end

%%
function [g_boundary_analyt, sig_scan]  = computeBoundary_productFeedback(beta1, beta2, p0, S10, S20)

% DOL steady state and reference biomass Md
sol_dol = solve_dol_productFeedback(beta1,beta2,p0,S10,S20);
Md = sol_dol(1,3) + sol_dol(1,4);

% Define phi, Phi1/2 and analytical gFun(sigma, Ms)
Phi1 = @(Ms) p0.phi(p0.ths.*Ms, beta1);
Phi2 = @(Ms) p0.phi(p0.ths.*Ms, beta2);

gFun = @(sig, Ms) p0.dilution /(p0.Ad1*p0.dilution + p0.md1) .* sig .*...
    (p0.k1.*S10.*Phi1(Ms)./(p0.dilution + Ms.*sig.*p0.k1.*Phi1(Ms)) + ...
     p0.k2.*S20.*Phi2(Ms)./(p0.dilution + Ms.*sig.*p0.k2.*Phi2(Ms)));

% Analytical equal-biomass boundary g(sigma; Ms = Md)

sig_scan = linspace(0,1,101);
g_boundary_analyt = gFun(sig_scan, Md);

end

%%
function beta1_crit = compute_beta1_crit_DOL1washout(p0,S10,S20)
% For beta2=0, solve the critical value of beta1 that just washes out DOL1

% DOL2-only steady state (beta2 = 0)
Ud2 = p0.Ad2 + p0.md2/p0.dilution;
Sd2 = Ud2 * p0.dilution / p0.k2;
Md2 = (S20 - Sd2) / Ud2;

P_DOL2 = p0.thd2 * Md2;  % product from DOL2 only

% exp(beta1*Hill_P) * k1 * S10 = Ad1*gamma + md1
Hill_P = p0.Hill(P_DOL2);
rhs = (p0.Ad1 * p0.dilution + p0.md1)./(p0.k1*S10);
beta1_crit = log(rhs) / Hill_P; 
end