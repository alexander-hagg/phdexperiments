% metric: rank diversity
nTests = 100;
n = 10;
pSize = 32;
clear metric;
    
for j=1:nTests
    
    X = 0.1*rand(n,2);
    Y = 0.1*rand(n,2);
    
    XFit = rand(n,1);    
    YFit = 2*rand(n,1);
    
    metric(j,:) = metricRRD(X, Y, XFit, YFit);
    %metric(j,:) = metricRRD(Y, X, YFit, XFit);
end

figure(2)
boxplot(metric);
grid on;
axis([0.5 2.5 0 1]);


