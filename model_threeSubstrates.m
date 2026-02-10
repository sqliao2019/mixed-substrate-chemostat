function dydt = model_threeSubstrates(~, y, p, S10, S20, S30)
%%%% ODE function

% Variables
Sd1 = y(1);
Sd2 = y(2);
Sd3 = y(3);
Ss1 = y(4);
Ss2 = y(5);
Ss3 = y(6);

if ~isa(y,'sym')
    Md1 = y(7) * (y(7)>p.cutoff);
    Md2 = y(8) * (y(8)>p.cutoff);
    Md3 = y(9) * (y(9)>p.cutoff);
    Ms  = y(10) * (y(10)>p.cutoff);
else
    Md1 = y(7);
    Md2 = y(8);
    Md3 = y(9);
    Ms  = y(10); 
end

% DOL
dSd1_dt = p.dilution * (S10 - Sd1) - p.f1(Sd1) * Md1;
dSd2_dt = p.dilution * (S20 - Sd2) - p.f2(Sd2) * Md2;
dSd3_dt = p.dilution * (S30 - Sd3) - p.f3(Sd3) * Md3;
Lamd1 = (p.f1(Sd1) - p.md1)/p.Ad1;
Lamd2 = (p.f2(Sd2) - p.md2)/p.Ad2;
Lamd3 = (p.f3(Sd3) - p.md3)/p.Ad3;
dMd1_dt = (Lamd1 - p.dilution) * Md1;
dMd2_dt = (Lamd2 - p.dilution) * Md2;
dMd3_dt = (Lamd3 - p.dilution) * Md3;

% SS
dSs1_dt = p.dilution * (S10 - Ss1) - p.sig1 * p.f1(Ss1) * Ms;
dSs2_dt = p.dilution * (S20 - Ss2) - p.sig2 * p.f2(Ss2) * Ms;
dSs3_dt = p.dilution * (S30 - Ss3) - p.sig3 * p.f3(Ss3) * Ms;
Lams = (p.sig1 * p.f1(Ss1) + p.sig2 * p.f2(Ss2) + p.sig3 * p.f3(Ss3) - p.ms) / p.As;
dMs_dt = (Lams - p.dilution) * Ms;

dydt = [dSd1_dt; dSd2_dt; dSd3_dt; ...
        dSs1_dt; dSs2_dt; dSs3_dt; ...
        dMd1_dt;  dMd2_dt;  dMd3_dt;  ...
        dMs_dt];

