tic

set(0,'DefaultAxesFontSize',20)

T = 2:2:800;
D = 0*T;

g = 9810;
gammaF = 4.29*g;
gamma = 1.047*gammaF;

for i = 1:length(T)
    D(i) = DiffCoeff(T(i),gamma);
end

Diff = median(D(end*3/4:end));
dDiff = std(D(end*3/4:end));

plot(T,D,'.','markersize',20)
hold on
plot([0 T(end)],[Diff Diff],'--k','linewidth',2)
hold on
gray = [0.8 0.8 0.8];
patch([0 T(end) T(end) 0],[Diff-dDiff Diff-dDiff Diff+dDiff Diff+dDiff],gray,'FaceAlpha',0.5)
%fill([0 T(end)], [Diff-dDiff Diff+dDiff],gray,'FaceAlpha',0.5)
hold off
xlabel('Elapsed time (t/T_f)')
ylabel('Approximate diffusivity (D)')
axis([0 T(end) 0 0.15])

toc