function y=gauss3(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3)

% gauss3
%
% y=gauss3(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3)
%
%Models a Gaussian curve to defined wavelength range using height, width and
%spectral location
%
% Inputs:
%     x        = wavelength values affiliated with Gaussian curves or fitting
%                spectral window
%     sigma1   = width of 1st Gaussian curve
%     phi1     = height of 1st Gaussian curve
%     mu1      = spectral location of 1st Gaussian curve (center)
%     sigma'n' = width of n-th Gaussian curve
%     phi'n'   = height of n-th Gaussian curve
%     mu'n'    = spectral location of n-th Gaussian curve (center)
%
% Returns:
%     y     = modeled Gaussian curve
% 
% copyright (c) 2018 Brice K. Grunert
% email: bricegrunert@gmail.com

%Last modified on 28 June 2018 by BG

y=phi1.*exp(-(x-mu1).^2./(2.*sigma1.^2))+phi2.*exp(-(x-mu2).^2./(2.*sigma2.^2))+phi3.*exp(-(x-mu3).^2./(2.*sigma3.^2));