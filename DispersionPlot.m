tic

set(0,'DefaultAxesFontSize',20)

x = [0:1/300:1 - 1/300];
T = 80;
g = 9810;
gammaF = 4.29*g;
gamma = 1.047*gammaF;

for i = 1:length(x)
    [taux] = KMmultICs_Dispersion(x(i), T, gamma);
    plot(taux(:,1),taux(:,2)-x(i))
    hold on
end

xlabel('Elapsed Time t/T_f')
ylabel('Displacement x/\lambda_f')
axis([0 80 -15 15])
hold off

toc