set(0,'DefaultAxesFontSize',20)

tic

taux = KM_Timeseries;

vx = taux(:, end-1);
vy = taux(:, end);
maxtau = max(taux(:, 3));

for i = 1:length(vx)
Col = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];
idx = round(taux(i,3)*6/maxtau + 1);
if idx < 1
idx = 1;
end 
plot(vx(i), vy(i), '.', 'Color', Col(idx, :), 'markersize', 0.5)
hold on
end

% h(1) = plot(NaN,NaN, '.', 'Color', Col(1, :), 'markersize', 40);
% h(2) = plot(NaN,NaN, '.', 'Color', Col(2, :), 'markersize', 40);
% h(3) = plot(NaN,NaN, '.', 'Color', Col(3, :), 'markersize', 40);
% h(4) = plot(NaN,NaN, '.', 'Color', Col(4, :), 'markersize', 40);
% h(5) = plot(NaN,NaN, '.', 'Color', Col(5, :), 'markersize', 40);
% h(6) = plot(NaN,NaN, '.', 'Color', Col(6, :), 'markersize', 40);
% h(7) = plot(NaN,NaN, '.', 'Color', Col(7, :), 'markersize', 40);
%legend(h, 'Quick bounce', '1/3 period', '2/3 period', '1 period', '4/3 periods', '5/3 periods', '2 periods','Location','EastOutside');

% xlabel('v_x in cm/s')
% ylabel('v_y in cm/s')
% axis([-20 20 -5 10])
xticks([])
yticks([])
ax.YColor = 'w'; % White
ax.XColor = 'w'; % White
axis off

box off

hold off

saveas(gcf, 'BengalTiger.fig')

toc