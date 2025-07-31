function g = ift2(G, delta_f)
% function g = ift2(G, delta_f)
% Numerical simulation of optical wave propagation with examples in Matlabs
% Jason D. Schmidt
% page 37
N = size(G, 1);
g = ifftshift(ifft2(ifftshift(G))) * (N * delta_f)^2;