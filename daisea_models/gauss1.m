function y=gauss1(x,sigma1,phi1,mu1)

% gauss1
%
% y=gauss1(x,sigma1,phi1,mu1,sigma2,phi2,mu2)
%
%Models a Gaussian curve to defined wavelength range using height, width and
%spectral location
%
% Inputs:
%     x        = wavelength values affiliated with Gaussian curves or fitting
%                spectral window
%     sigma1   = width of Gaussian curve
%     phi1     = height of Gaussian curve
%     mu1      = spectral location of Gaussian curve (center)
%
% Returns:
%     y     = modeled Gaussian curve
% 

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

%Last modified on 28 June 2018 by BG

y=phi1.*exp(-(x-mu1).^2./(2.*sigma1.^2));