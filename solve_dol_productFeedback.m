function sol_dol = solve_dol_productFeedback(beta1, beta2, p0, S10, S20)
%SOL_DOL_PRODUCTFEEDBACK  Steady state of the DOL subsystem with product feedback.
%
%   This function does NOT solve the full 5-equation DOL ODE subsystem
%   directly.  Instead, it uses analytically reduced steady-state relations
%   that are exact for the Fig. 4 parameter regimes:
%
%   (1) Symmetric case:
%       beta1 = beta2, S10 = S20, and
%       Ad1 = Ad2, md1 = md2, thd1 = thd2, k1 = k2.
%       Under these conditions, the symmetric DOL steady state
%           Sd1 = Sd2,  Md1 = Md2
%       can be determined by solving a single scalar equation for P, with
%       Sd and Md then obtained directly from algebraic formulas.
%
%   (2) Asymmetric case with beta2 = 0:
%       DOL2 (no product feedback) has an explicit closed-form steady state
%       obtained from its algebraic balance equations.
%       Given this baseline production from DOL2, the remaining DOL1+Pd
%       steady state is again found by solving a single scalar equation for
%       P, followed by algebraic back-substitution for Sd1 and Md1.
%
%   In both regimes, the only numerical root-finding step is a scalar solve
%   in P; all substrate and biomass values are obtained by direct formulas.
%
%   The function returns the DOL steady state in the order:
%       [Sd1  Sd2  Md1  Md2  Pd]

    tol = 1e-8;

    isSymBeta   = abs(beta1 - beta2) < tol;
    isSymInput  = abs(S10 - S20)    < tol;
    isSymParams = ...
        abs(p0.Ad1 - p0.Ad2)   < tol && ...
        abs(p0.md1 - p0.md2)   < tol && ...
        abs(p0.thd1 - p0.thd2) < tol && ...
        abs(p0.k1  - p0.k2)    < tol;

    % ---- Branch 1: Fully symmetric DOL ---- %
    if isSymBeta && isSymInput && isSymParams
        sol_dol = solve_dol_symmetric(beta1, p0, S10);
        return;
    end

    % ---- Branch 2: beta2 = 0 (DOL2 independent of product) ---- %
    if abs(beta2) < tol
        sol_dol = solve_dol_beta2_zero(beta1, p0, S10, S20);
        return;
    end

    error('solve_dol_productFeedback is only implemented for beta1=beta2 or beta2=0.');
end


% =================== Symmetric case beta1 = beta2 =================== %

function sol_dol = solve_dol_symmetric(beta1, p0, S10)
% SOLVE_DOL_SYMMETRIC
%   Steady state of the symmetric DOL subsystem (beta1 = beta2).
%   Uses growth = dilution instead of dM/dt = 0, eliminating the trivial branch.
    tolS = 1e-10;
    gamma = p0.dilution;
    Ud1   = p0.Ad1 + p0.md1/gamma;
    
    if Ud1 <= 0
        error('Metabolic cost must be positive');
    end

    if S10 < -tolS
        error('Substrate input must be non-negative');
    end

    % search interval for nontrivial steady-state
    Md1_min = 0;
    Md1_max = S10 / Ud1;  

    if Md1_max <= 0
        % S10 = 0
        Md1_sol = 0;
        Sd1_sol = S10;
        Pd_sol = 0;
    else
        syms Md1 real   % biomass of each DOL strain (identical)
    
        % product balance Pd = thd1*Md1 + thd2*Md2 = 2*thd1*Md1 (symmetric)
        Pd_from_Md1  = @(M) 2 * p0.thd1 * M;
    
        % substrate balance Sd1 = S10 - Ud1 * Md1
        Sd1_from_Md1 = @(M) S10 - Ud1 * M;
    
        % growth = dilution  → phi(Pd) * f1(Sd1) = Ud1 * gamma
        eqn = Ud1 * gamma == ...
              p0.phi(Pd_from_Md1(Md1), beta1) .* p0.f1(Sd1_from_Md1(Md1));
    
        Md1_sol = vpasolve(eqn, Md1, [Md1_min Md1_max]);
    
        if isempty(Md1_sol)
            Md1_sol = 0;   % trivial washout
        else
            Md1_sol = double(Md1_sol);
        end
    
        Sd1_sol = Sd1_from_Md1(Md1_sol);
        Pd_sol  = Pd_from_Md1(Md1_sol);
    end
    % symmetric duplication
    sol_dol = [Sd1_sol, Sd1_sol, Md1_sol, Md1_sol, Pd_sol];
end


% ==================== Asymmetric case beta2 = 0 ===================== %

function sol_dol = solve_dol_beta2_zero(beta1, p0, S10, S20)
    Stot = S10 + S20; 
    tolS = 1e-10 * max(1, abs(Stot));
    
    gamma = p0.dilution;
    Ud1 = p0.Ad1 + p0.md1/gamma;
    Ud2 = p0.Ad2 + p0.md2/gamma;
    if Ud1 <= 0 || Ud2 <= 0
        error('Metabolic cost must be positive');
    end   
    if S10 < - tolS || S20 < - tolS
        error('Substrate input must be non-negative');
    end
    % ---- DOL2: explicit solution (beta2 = 0) ----
    % search interval for nontrivial steady-state
    Md2_min = 0;
    Md2_max = S20 / Ud2;
    
    if Md2_max <= 0
        % S20 = 0
        Md2_sol = 0;
        Sd2_sol = S20;
    else
        syms Md2 real
        Sd2_from_Md2 = @(M) S20 - Ud2*M;
        
        % For beta2 = 0 and phi(P,0) = 1, growth = dilution reduces to f2(Sd2) = Ud2*gamma
        eqn_Md2 = Ud2 * gamma == p0.f2(Sd2_from_Md2(Md2));
    
        Md2_sol = vpasolve(eqn_Md2, Md2, [Md2_min Md2_max]);
    
        if isempty(Md2_sol)
            Md2_sol = 0;   % trivial washout
        else
            Md2_sol = double(Md2_sol);
        end
        Sd2_sol = Sd2_from_Md2(Md2_sol);
    end

    baseP = p0.thd2 * Md2_sol;  % baseline product from DOL2

    % ---- DOL1 + Pd: scalar equation in Md1 ----
    % search interval for nontrivial steady-state
    Md1_min = 0;
    Md1_max = S10 / Ud1;

    if Md1_max <= 0
        % S10 = 0
        Md1_sol = 0;
        Sd1_sol = S10;
        Pd_sol = baseP;
    else
        syms Md1 real   % biomass of DOL1
    
        % product balance Pd = thd1*Md1 + thd2*Md2
        Pd_from_Md1  = @(M) p0.thd1 * M + baseP;
    
        % substrate balance Sd1 = S10 - Ud1 * Md1
        Sd1_from_Md1 = @(M) S10 - Ud1 * M;
    
        % growth = dilution  → phi(Pd) * f1(Sd1) = Ud1 * gamma
        eqn = Ud1 * gamma == ...
              p0.phi(Pd_from_Md1(Md1), beta1) .* p0.f1(Sd1_from_Md1(Md1));
    
    
        Md1_sol = vpasolve(eqn, Md1, [Md1_min Md1_max]);
    
        if isempty(Md1_sol)
            Md1_sol = 0;   % trivial washout
        else
            Md1_sol = double(Md1_sol);
        end
    
        Sd1_sol = Sd1_from_Md1(Md1_sol);
        Pd_sol  = Pd_from_Md1(Md1_sol);
    end
    sol_dol = [Sd1_sol, Sd2_sol, Md1_sol, Md2_sol, Pd_sol];
end
