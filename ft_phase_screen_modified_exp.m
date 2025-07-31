function phz = ft_phase_screen_modified(r0, N, delta, L0, l0)
% Numerical simulation of optical wave propagation with examples in Matlabs
% Jason D. Schmidt
% page 168

% function phz ...
% = ft_phase_screen(r0, N, delta, L0, l0)

% setup the PSD
del_f = 1/(N*delta); % frequency grid spacing [1/m]
fx = (-N/2 : N/2-1) * del_f;
% frequency grid [1/m]
[fx,fy] = meshgrid(fx);
[th,f] = cart2pol(fx, fy); % polar grid
fl = 3.3/l0/(2*pi); % inner scale frequency [1/m]
f0 = 4/L0; % outer scale frequency [1/m]
% modified von Karman atmospheric phase PSD
PSD_phi = 0.023*r0^(-5/3) * (1+1.802*(f/fl)-0.254*(f/fl).^(7/6)).*(1-exp(-f.^2/f0.^2)).*exp(-f.^2/fl.^2)./f.^(11/3);%/((f.^2+f0.^2).^(11/6));
PSD_phi(N/2+1,N/2+1) = 0;
% random draws of Fourier coefficients
cn = (randn(N) + 1i*randn(N)) .* sqrt(PSD_phi)*del_f;
% synthesize the phase screen
phz = real(ift2(cn, 1)); % 
% % % phz = real(ifft2(cn, 1)); % ? why sppecify row number 1
% % % phz = real(ifft2(cn)); % has aliasing probCDY5    tc Cqlem, pi jump
% phz = real(ifftshift(ifft2(ifftshift(cn))));