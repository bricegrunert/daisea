function y=cdom_model_exponential(x,alam0,s,lam0)

% y=cdom_model_exponential(x,alam0,s,lam0)
%
% Models an exponential function to fit CDOM and/or NAP absorption
% spectra using wavelength, absorption at an initial wavelength and 
% spectral slope
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

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

%Last modified on 28 June 2018 by BG


y=alam0.*exp(-(x-lam0).*s);