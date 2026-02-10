function dydt = model_general(~, y, p, S10, S20)
%%%% ODE function

% Variables
Sd1 = y(1);
Sd2 = y(2);
Ss1 = y(3);
Ss2 = y(4);

if ~isa(y,'sym')
    Md1 = y(5) * (y(5)>p.cutoff);
    Md2 = y(6) * (y(6)>p.cutoff);
    Ms = y(7) * (y(7)>p.cutoff);
else
    Md1 = y(5);
    Md2 = y(6);
    Ms = y(7); 
end

% DOL
dSd1_dt = p.dilution * (S10 - Sd1) - p.f1(Sd1) * Md1;
dSd2_dt = p.dilution * (S20 - Sd2) - p.f2(Sd2) * Md2;
Lamd1 = (p.f1(Sd1) - p.md1)/p.Ad1;
Lamd2 = (p.f2(Sd2) - p.md2)/p.Ad2;
dM1_dt = (Lamd1 - p.dilution) * Md1;
dM2_dt = (Lamd2 - p.dilution) * Md2;

% SS
dS1s_dt = p.dilution * (S10 - Ss1) - p.sig1 * p.f1(Ss1) * Ms;
dS2s_dt = p.dilution * (S20 - Ss2) - p.sig2 * p.f2(Ss2) * Ms;
Lams = (p.sig1 * p.f1(Ss1) + p.sig2 * p.f2(Ss2) - p.ms) / p.As;
dMs_dt = (Lams - p.dilution) * Ms;

dydt = [dSd1_dt; dSd2_dt; dS1s_dt; dS2s_dt; dM1_dt; dM2_dt; dMs_dt];