function y=cdom_model_exponential_k(x,a,s,lam0,k)

% y=cdom_model_exponential_k(x,a,s,lam0,k)
%
% Exponential model for CDOM slope.
%
% Inputs:
%     x     = wavelength values affiliated with absorption spectra
%     a     = CDOM and/or NAP absorption at lam0
%     s     = CDOM and/or NAP absorptio spectral slope
%     lam0  = wavelength used to initialize fitting
%     k     = offset
%
% Returns:
%     y     = modeled absorption spectra
% 


% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu


y=a.*exp(-(x-lam0).*s)+k;