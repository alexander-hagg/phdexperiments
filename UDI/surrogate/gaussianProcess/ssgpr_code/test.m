clear;clc;
numSpectralPts = 50;
numTrSamples = 2000; numDim = 1;
numTestSamples = 10000;
fcn = @(X_tr) sin(X_tr) + sin(2*X_tr).^2 + 5*(sin(10*X_tr).^6.*tanh(0.02*X_tr));
X_tr    = 2*pi*rand(numTrSamples,numDim);
X_tr    = sort(X_tr);
T_tr    = fcn(X_tr);
X_tst   = 2*pi*rand(numTestSamples,numDim);
X_tst   = sort(X_tst);
T_tst   = fcn(X_tst);

%
[NMSE, mu, S2, NMLP, loghyper, convergence, lengthscales] = ssgpr_ui(X_tr, T_tr, X_tst, T_tst, numSpectralPts); 
%
figure(1);hold off;
plot(X_tr,T_tr,'x');
hold on;
plot(X_tst,T_tst,'-');
plot(X_tst,mu,'-');
fill([X_tst;flipud(X_tst)],[mu+S2;flipud(mu-S2)],'-');
alpha(0.25);
drawnow;
legend('Samples', 'Target', 'Prediction');
