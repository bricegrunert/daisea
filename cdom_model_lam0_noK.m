function y=cdom_model_lam0_noK(x,alam0,s,lam0)

% cdom_model_lam0_noK
%
% y=cdom_model_lam0_noK(x,a,s,lam0)
%
%Models an exponential function to fit CDOM and/or NAP absorption
%spectra using wavelength, absorption at an initial wavelength and spectral
%slope
%
% Inputs:
%     x     = wavelength values affiliated with absorption spectra
%     lam0  = wavelength used to initialize fitting
%     alam0 = CDOM and/or NAP absorption at lam0
%     s     = CDOM and/or NAP absorptio spectral slope
%
% Returns:
%     y     = modeled absorption spectra
% 
% copyright (c) 2018 Brice K. Grunert
% email: bricegrunert@gmail.com

%Last modified on 28 June 2018 by BG


y=alam0.*exp(-(x-lam0).*s);