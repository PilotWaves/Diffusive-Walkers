tic

lam_f = 4.75; % dim wavelength

dx = 0.01;
x = [0:dx:lam_f - dx];
DX = 0*x;
N = 10000;
XX = zeros(N,length(x));

parfor i = 1:length(x)
    [X] = KMmultICs(x(i), N)
    DX(i) = X(end) - x(i);
    XX(:,i) = (X - x(i))/lam_f;
end

plot([1:N],XX)

D = mean(DX.^2/lam_f^2)/2/N %change N to time

toc