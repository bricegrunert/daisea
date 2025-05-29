function y=gauss(x,sigma,phi,mu)

% gauss
%
% y=gauss(x,sigma,phi,mu)
%
%Models a Gaussian curve to defined wavelength range using height, width and
%spectral location
%
% Inputs:
%     x     = wavelength values affiliated with Gaussian curve
%     sigma = width of Gaussian curve
%     phi   = height of Gaussian curve
%     mu    = spectral location of Gaussian curve (center)
%
% Returns:
%     y     = modeled Gaussian curve
% 

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

%Last modified on 28 June 2018 by BG


y=phi.*exp(-(x-mu).^2./(2.*sigma.^2));