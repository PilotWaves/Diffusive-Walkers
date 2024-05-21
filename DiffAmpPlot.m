tic

set(0,'DefaultAxesFontSize',20)

T = 80;

g = 9810;
gammaF = 4.29*g;
%sqrteps = [0, 0.04, 0.07, 0.08, 0.09, 0.12, 0.16, 0.18, 0.2, sqrt(0.047)];
gamma = [1:.001:1.047]*gammaF;
D = 0*gamma;

for i = 1:length(gamma)
    D(i) = DiffCoeff(T, gamma(i));
end

%eps = gamma/gammaF - 1; % reduced acceleration

% a = 0.49/sqrt(1.047 - 1); % dimensional amplitude at gamma = 1.047gammaF
% A = a*sqrt(eps); % dimensional amplitude

plot(gamma/gammaF,D,'.','markersize',20)
% hold on
% plot(sqrteps,D(1:length(sqrteps)),'s','markersize',40)
% hold off
xlabel('Root reduced acceleration $(\sqrt{\gamma/\gamma_F - 1})$','Interpreter','latex')
ylabel('Approximate diffusivity (D)')
axis([0 0.23 0 0.15])

toc