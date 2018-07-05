function [model, model_gof, model_fitstats]=spectral_slope_nogauss_noK(absorption,wavelength,lambda_start,lambda_stop,lam0)

% spectral_slope_nogauss_noK
%
% [model, model_gof, model_fitstats]=spectral_slope_nogauss_noK(absorption,wavelength,lambda_start,lambda_stop,lam0)
%
%Calculates spectral slope for CDOM data within a defined spectral range
%Ranges can go from 350-500nm (standard spectral slope; Stedmon et al.
%2000), 240-700nm + Gaussian components (Massicotte & Markager 2016), 
%275-295nm and 350-400nm and their ratio (Sr,dimensionless) to provide 
%insight on source (marine or terrestrial; Helms et al. 2008)
%
% Inputs:
%     absorption   = cdom/nap absorption spectra, format = vector
%     wavelength   = wavelength values affiliated with absorption spectra,
%                    format = vector
%     lambda_start = minimum wavelength considered (e.g. 350 nm),
%                    must be within spectral range of wavelength input
%     lambda_stop  = maximum wavelength considered (e.g. 700 nm),
%                    must be within spectral range of wavelength input
%     lam0         = wavelength used to initialize fitting
%
% Returns:
%     model          = output from fitting, including spectral slope
%                      as model.s and absorption at lam0 as model.a
%     model_gof      = goodness of fit for model
%     model_fitstats = fitting stats for model
% 
% copyright (c) 2018 Brice K. Grunert
% email: bricegrunert@gmail.com

%Last modified on 28 June 2018 by BG

%%

%Index matching input wavelength to specified minimum and maximum
%wavelength
for jj=1:length(wavelength)
    if wavelength(jj)==lambda_start
        ind=jj;
    end
    if wavelength(jj)==lambda_stop
        ind2=jj;
    end
end

if size(absorption,2)>1
    absorption=absorption';
end

if size(wavelength,2)>1
    wavelength=wavelength';
end

x=wavelength(ind:ind2);
y=absorption(ind:ind2);

ind=find(isnan(y)==0);

x=x(ind);
y=y(ind);

ft=fittype('cdom_model_lam0_noK(x,a,s,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s'},'problem',{'lam0'}); % to make something a constant, specify it as 'problem',{'variable'}
fopts=fitoptions(ft);
fopts.StartPoint=[5 0.01];

[model,model_gof,model_fitstats]=fit(x,y,ft,fopts,'problem',{lam0}); % define problem variables, if any at the end as 'problem',{value1,value2}
