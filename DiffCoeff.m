function [D] = DiffCoeff(T, gamma)

x = [0:1/300:1 - 1/300];
DX = 0*x;

parfor i = 1:length(x)
    [X] = KMmultICs(x(i), T, gamma);
    DX(i) = X - x(i);
end

    D = mean(DX.^2)/(2*T); %replace N with time