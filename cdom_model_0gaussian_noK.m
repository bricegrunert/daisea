function y=cdom_model_0gaussian_noK(x,a,s,lam0)

% cdom_model_0gaussian_noK
%
% y=cdom_model_0gaussian_noK(x,a,s1,lam0)
%
% Models a single exponential with no Gaussian curves defined
%
% Inputs:
%     x        = wavelength values affiliated with Gaussian curves or fitting
%                spectral window
%     lam0     = wavelength used to initialize fitting
%     alam0    = CDOM and/or NAP absorption at lam0
%     s        = CDOM and/or NAP absorptio spectral slope
%
% Returns:
%     y     = modeled spectra
% 
% copyright (c) 2018 Brice K. Grunert
% email: bricegrunert@gmail.com

%Last modified on 28 June 2018 by BG

y=a.*exp(-(x-lam0).*s);