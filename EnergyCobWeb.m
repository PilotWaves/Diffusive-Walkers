set(0,'DefaultAxesFontSize',20)

rho = 9.5e-7; % dim oil density
R = 0.28; % dim drop radius
m = rho*R^3;
u_in = m*taux(:, end-1).^2;
u_out = m*taux(:, end).^2;

tic
for i = 1:length(u_in)-2
% Col = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];
% idx = floor(log10(abs(v(i+1)-v(i)))+length(Col)-1);
% if idx < 1
% idx = 1;
% end
%plot([abs(v(i)) abs(v(i))], [abs(v(i)) abs(v(i+1))], [abs(v(i)) abs(v(i+1))], [abs(v(i+1)) abs(v(i+1))], '--', 'Color', Col(idx, :), 'linewidth', 0.5)
plot([abs(u_in(i)) abs(u_out(i))], [abs(u_in(i+1)) abs(u_out(i+1))], [abs(u_out(i)) abs(u_in(i+1))], [abs(u_out(i+1)) abs(u_in(i+2))], 'linewidth', 0.5)
hold on
end
xlabel('v_{in} in cm/s')
ylabel('v_{out} in cm/s')
hold off
toc