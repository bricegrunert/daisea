function y=cdom_model_stretched_exponential(x,a,s,beta,lam0)

% y=cdom_model_stretched_exponential(x,a,s,beta,lam0)
%
% Generic stretched exponential, based off Cael & Boss 2017
% beta is a value from 0-1
%
% Inputs:
%     x     = wavelength values affiliated with absorption spectra
%     a     = CDOM and/or NAP absorption at lam0
%     s     = CDOM and/or NAP absorptio spectral slope
%     beta  = stretch parameter
%     lam0  = wavelength used to initialize fitting
%
% Returns:
%     y     = modeled absorption spectra
% 

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu


y=a.*exp(-(s.*(x-lam0)).^beta);