function y=cdom_model_hyperbolic(x,a,s)

% y=cdom_model_hyperbolic(x,a,s)
% 
% Hyperbolic model for fitting CDOM/NAP. Twardowski et al. 2004 Marine
% Chemistry.
%
% Inputs:
%     x     = wavelength values affiliated with absorption spectra
%     a     = CDOM and/or NAP absorption at 440 nm
%     s     = CDOM and/or NAP absorptio spectral slope
%
% Returns:
%     y     = modeled absorption spectra
% 


% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

% Created on 2022-12-01 by AC

%%
%y=a.*(x/532).^-s;
y=a.*(x/440).^-s;