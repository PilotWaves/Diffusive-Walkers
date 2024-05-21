set(0,'DefaultAxesFontSize',20)

v = taux(:, end-1);

tic
for i = 1:length(v)-1
Col = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];
idx = round(taux(i,3)*6/maxtau + 1);
if idx < 1
idx = 1;
end 
plot([abs(v(i)) abs(v(i))], [abs(v(i)) abs(v(i+1))], [abs(v(i)) abs(v(i+1))], [abs(v(i+1)) abs(v(i+1))], '--', 'Color', Col(idx, :), 'linewidth', 0.5)
%plot([abs(v(i)) abs(v(i))], [abs(v(i)) abs(v(i+1))], [abs(v(i)) abs(v(i+1))], [abs(v(i+1)) abs(v(i+1))], 'linewidth', 0.5)
hold on
end

h(1) = plot(NaN,NaN, 'Color', Col(1, :), 'linewidth', 10);
h(2) = plot(NaN,NaN, 'Color', Col(2, :), 'linewidth', 10);
h(3) = plot(NaN,NaN, 'Color', Col(3, :), 'linewidth', 10);
h(4) = plot(NaN,NaN, 'Color', Col(4, :), 'linewidth', 10);
h(5) = plot(NaN,NaN, 'Color', Col(5, :), 'linewidth', 10);
h(6) = plot(NaN,NaN, 'Color', Col(6, :), 'linewidth', 10);
h(7) = plot(NaN,NaN, 'Color', Col(7, :), 'linewidth', 10);
legend(h, 'Quick bounce', '1/3 period', '2/3 period', '1 period', '4/3 periods', '5/3 periods', '2 periods','Location','EastOutside');


xlabel('v_{in} in cm/s')
ylabel('v_{out} in cm/s')
axis([0 20 0 20])
hold off
toc