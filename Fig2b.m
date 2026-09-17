Param = getPara(5);

numVF = 16;
x1v = linspace(0, Param.S1in/Param.a1, numVF);
x2v = linspace(0, Param.S2in/Param.a2, numVF);

[X1,X2] = meshgrid(x1v, x2v);

S1 = Param.S1in - Param.a1 .* X1;
S2 = Param.S2in - Param.a2 .* X2;
M1 = Param.M1in + Param.k1 .* X1 - Param.b2 .* X2;
M2 = Param.M2in + Param.k2 .* X2 - Param.b1 .* X1;

% Numerical tolerance for feasibility checks to avoid
% floating-point exclusion of points on the boundary
tol = 1e-8;
mask = S1<-tol | S2<-tol | M1<-tol | M2<-tol;

Fx1 = (Param.U1(S1,M2) - Param.dilution).* X1;
Fx2 = (Param.U2(S2,M1) - Param.dilution).* X2;

Fx1(mask) = NaN;
Fx2(mask) = NaN;

mag = sqrt(Fx1.^2 + Fx2.^2);

Fx1_norm = Fx1./mag;
Fx2_norm = Fx2./mag;

magNorm = (mag - min(mag(:))) ./ (max(mag(:)) - min(mag(:)));
cmap = parula(256);

figure;
axis equal; hold on;

for i = 1:numel(X1)
    % Pick color based on magnitude
    idx = max(1, round(magNorm(i)*255)+1);
    thisColor = cmap(idx,:);
    
    % Plot each arrow with that color
    h = quiver(X1(i), X2(i), Fx1_norm(i), Fx2_norm(i), 0.15, ...
           'Color', thisColor, 'MaxHeadSize', 3,'LineWidth', 1.2);
    h.ShowArrowHead = 'on'; 

end

edgewidth = 0.00;
axis([0-edgewidth, Param.S1in/Param.a1+edgewidth, 0-edgewidth, Param.S2in/Param.a2+edgewidth]);

colormap(cmap);
cb = colorbar;
cb.Label.String = 'Vector magnitude';

plot(x1v, (Param.b1 .*x1v - Param.M2in) ./Param.k2 , 'k');
plot((Param.b2 .*x2v - Param.M1in) ./Param.k1 , x2v, 'k');


% Fixed points and separatrices
fixedPoint = solveFixedPoint(Param);
for i = 1:size(fixedPoint,1)

    xy = fixedPoint{i,1};

    % Plot steady state according to stability
    switch fixedPoint{i,2}
        case 1          % stable
            plot(xy(1), xy(2), 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.1 0.1 0.1]);
        case 2          % non-hyperbolic
            plot(xy(1), xy(2), 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[0.5 0.5 0.5]);
        case 0          % unstable
            plot(xy(1), xy(2), 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 1 1]);
            % Plot separatrix for saddle point
            if ~isempty(fixedPoint{i,3})
                plot(fixedPoint{i,3}(:,1), fixedPoint{i,3}(:,2), 'k');
            end
    end


end

Data = [X1(:),X2(:),Fx1(:),Fx2(:),mag(:),Fx1_norm(:),Fx2_norm(:)];