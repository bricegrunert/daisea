function y=cdom_model_3gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,lam0)

% cdom_model_3gaussian_noK
%
% y=cdom_model_3gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,lam0)
%
% Models a single exponential and 3 Gaussian curves to estimate an
% absorption spectra
%
% Inputs:
%     x        = wavelength values affiliated with Gaussian curves or fitting
%                spectral window
%     lam0     = wavelength used to initialize fitting
%     alam0    = CDOM and/or NAP absorption at lam0
%     s        = CDOM and/or NAP absorptio spectral slope
%     sigma1   = width of 1st Gaussian curve
%     phi1     = height of 1st Gaussian curve
%     mu1      = spectral location of 1st Gaussian curve (center)
%     sigma'n' = width of n-th Gaussian curve
%     phi'n'   = height of n-th Gaussian curve
%     mu'n'    = spectral location of n-th Gaussian curve (center)
%
% Returns:
%     y     = modeled spectra
% 
% copyright (c) 2018 Brice K. Grunert
% email: bricegrunert@gmail.com

%Last modified on 28 June 2018 by BG

y=a.*exp(-(x-lam0).*s)+phi1.*exp(-(x-mu1).^2./(2.*sigma1.^2))+phi2.*exp(-(x-mu2).^2./(2.*sigma2.^2))+phi3.*exp(-(x-mu3).^2./(2.*sigma3.^2));