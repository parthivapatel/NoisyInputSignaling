Agonist = [.005,.01,.05,.1,.5,1];
atRA = [2/103,17/96,44/119,46/95,.72,.82,];

% MAPPING: Emax = b(1),  EC50 = b(2)
hill_fit = @(b,x)  x.^b(2)./(b(1)+x.^b(2));

b0 = [10,.5];                                  % Initial Parameter Estimates

B = lsqcurvefit(hill_fit, b0, Agonist, atRA);

AgVct = linspace(10^-3, 10,10000);   % Plot Finer Resolution

figure(1)
plot(Agonist, atRA, 'o')
hold on
plot(AgVct, hill_fit(B,AgVct), '-r')
hold off
legend('Data', 'Hill Equation Fit (n=.7)', 'Location','SE')
