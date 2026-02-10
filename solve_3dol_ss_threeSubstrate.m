function sol = solve_3dol_ss_threeSubstrate(p0, S10, S20, S30)
    
    gamma = p0.dilution;
    Stot = S10 + S20 + S30; 
    tolS = 1e-10 * max(1, abs(Stot));

%% ---- DOL1: explicit solution ----
    Ud1 = p0.Ad1 + p0.md1/gamma;
    if Ud1 <=0
        error('Metabolic cost must be positive');
    end
    if S10 < - tolS
        error('Substrate input must be non-negative');
    end

    Md1_min = 0;
    Md1_max = S10 / Ud1;  
    
    if Md1_max <= 0
        % S10 = 0
        Md1_sol = 0;
        Sd1_sol = S10;        
    else       
        syms Md1 real
        Sd1_from_Md1 = @(M) S10 - Ud1*M;
        
        % f1(Sd1) = Ud1*gamma
        eqn_Md1 = Ud1 * gamma == p0.f1(Sd1_from_Md1(Md1));
    
        Md1_sol = vpasolve(eqn_Md1, Md1, [Md1_min Md1_max]);
    
        if isempty(Md1_sol)
            Md1_sol = 0;   % trivial washout
        else
            Md1_sol = double(Md1_sol);
        end
        Sd1_sol = Sd1_from_Md1(Md1_sol);
    end

    
%% ---- DOL2: explicit solution ----
    Ud2 = p0.Ad2 + p0.md2/gamma;
    if Ud2 <= 0
        error('Metabolic cost must be positive');
    end
    if S20 < - tolS
        error('Substrate input must be non-negative');
    end
    
    Md2_min = 0;
    Md2_max = S20 / Ud2;
    
    if Md2_max <= 0
        % S20 = 0: washout
        Md2_sol = 0;
        Sd2_sol = S20;
    else
        syms Md2 real
        Sd2_from_Md2 = @(M) S20 - Ud2*M;
        
        % f2(Sd2) = Ud2*gamma
        eqn_Md2 = Ud2 * gamma == p0.f2(Sd2_from_Md2(Md2));
        
        Md2_sol = vpasolve(eqn_Md2, Md2, [Md2_min Md2_max]);
        if isempty(Md2_sol)
            Md2_sol = 0;   % trivial washout
            Sd2_sol = S20;
        else
            Md2_sol = double(Md2_sol);
            Sd2_sol = Sd2_from_Md2(Md2_sol);
        end
    end

%% ---- DOL3: explicit solution ----
    Ud3 = p0.Ad3 + p0.md3/gamma;
    if Ud3 <= 0
        error('Metabolic cost must be positive');
    end
    if S30 < - tolS
        error('Substrate input must be non-negative');
    end
    
    Md3_min = 0;
    Md3_max = S30 / Ud3;
    
    if Md3_max <= 0
        % S30 = 0: washout
        Md3_sol = 0;
        Sd3_sol = S30;
    else
        syms Md3 real
        Sd3_from_Md3 = @(M) S30 - Ud3*M;
        
        % f3(Sd3) = Ud3*gamma
        eqn_Md3 = Ud3 * gamma == p0.f3(Sd3_from_Md3(Md3));
        
        Md3_sol = vpasolve(eqn_Md3, Md3, [Md3_min Md3_max]);
        if isempty(Md3_sol)
            Md3_sol = 0;   % trivial washout
            Sd3_sol = S30;
        else
            Md3_sol = double(Md3_sol);
            Sd3_sol = Sd3_from_Md3(Md3_sol);
        end
    end

%% ----- SS: explicit solution -----
    Us = p0.As + p0.ms/gamma;
    if Us <= 0
        error('Metabolic cost must be positive');
    end
    if S10 < - tolS || S20 < - tolS || S30 < - tolS
        error('Substrate input must be non-negative');
    end
    
    Ms_min = 0;
    Ms_max = (S10 + S20 + S30)/Us;

    if Ms_max<=0
        % no input
        Ms_sol = 0;
        Ss1_sol = S10;
        Ss2_sol = S20;
        Ss3_sol = S30;
    else
        syms Ms real
        Ss1 = SS_substrate_from_Ms(Ms, S10, p0.sig1, p0.alpha1, p0.K1, gamma);
        Ss2 = SS_substrate_from_Ms(Ms, S20, p0.sig2, p0.alpha2, p0.K2, gamma);
        Ss3 = SS_substrate_from_Ms(Ms, S30, p0.sig3, p0.alpha3, p0.K3, gamma);
    
        % sig1*f1(Ss1)+sig2*f2(Ss2)+sig3*f3(Ss3) = Us*gamma
        eqn_Ms = Us*gamma == p0.sig1*p0.f1(Ss1) + p0.sig2*p0.f2(Ss2) + p0.sig3*p0.f3(Ss3);
    
        Ms_sol = vpasolve(eqn_Ms, Ms, [Ms_min Ms_max]);
    
        if isempty(Ms_sol)
            Ms_sol = 0;   % trivial washout
        else
            Ms_sol = double(Ms_sol);
        end    
        
        Ss1_sol = SS_substrate_from_Ms(Ms_sol, S10, p0.sig1, p0.alpha1, p0.K1, gamma);
        Ss2_sol = SS_substrate_from_Ms(Ms_sol, S20, p0.sig2, p0.alpha2, p0.K2, gamma);
        Ss3_sol = SS_substrate_from_Ms(Ms_sol, S30, p0.sig3, p0.alpha3, p0.K3, gamma);
    end

%% ------------ output ------------
    sol = [Sd1_sol, Sd2_sol, Sd3_sol, ...
           Ss1_sol, Ss2_sol, Ss3_sol, ...
           Md1_sol, Md2_sol, Md3_sol, Ms_sol];
end

function S = SS_substrate_from_Ms(Ms, S0, sig, alpha, K, gamma)

    a = gamma;
    b = sig*alpha*Ms - gamma*(S0 - K);
    c = -gamma*S0*K;

    disc = b.^2 - 4*a*c;

    % Positive root
    S = (-b + sqrt(disc)) / (2*a);
end