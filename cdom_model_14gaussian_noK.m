function y=cdom_model_14gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13,sigma14,phi14,mu14,lam0)

% cdom_model_14gaussian_noK
%
% y=cdom_model_14gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13,sigma14,phi14,mu14,lam0)
%
% Models a single exponential and 14 Gaussian curves to estimate an
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

y=a.*exp(-(x-lam0).*s)+phi1.*exp(-(x-mu1).^2./(2.*sigma1.^2))+phi2.*exp(-(x-mu2).^2./(2.*sigma2.^2))+phi3.*exp(-(x-mu3).^2./(2.*sigma3.^2))+phi4.*exp(-(x-mu4).^2./(2.*sigma4.^2))+phi5.*exp(-(x-mu5).^2./(2.*sigma5.^2))+phi6.*exp(-(x-mu6).^2./(2.*sigma6.^2))+phi7.*exp(-(x-mu7).^2./(2.*sigma7.^2))+phi8.*exp(-(x-mu8).^2./(2.*sigma8.^2))+phi9.*exp(-(x-mu9).^2./(2.*sigma9.^2))+phi10.*exp(-(x-mu10).^2./(2.*sigma10.^2))+phi11.*exp(-(x-mu11).^2./(2.*sigma11.^2))+phi12.*exp(-(x-mu12).^2./(2.*sigma12.^2))+phi13.*exp(-(x-mu13).^2./(2.*sigma13.^2))+phi14.*exp(-(x-mu14).^2./(2.*sigma14.^2));