function [model,model_gof,model_fitstats]=daisea_fit_cdom_nap(absorption,wavelength,cdom_model,lambda_start,lambda_stop,lam0)

% [model,model_gof,model_fitstats]=daisea_fit_cdom_nap(absorption,wavelength,cdom_model,lambda_start,lambda_stop,lam0)
%
% Fits CDOM/NAP models for data within a defined spectral range. Ranges can
% go from 350-500nm (standard spectral slope; Stedmon et al. 2000), 
% 240-700nm + Gaussian components (Massicotte & Markager 2016), 
% 275-295nm and 350-400nm and their ratio (Sr,dimensionless) to provide 
% insight on source (marine or terrestrial; Helms et al. 2008)
%
% Inputs:
%     absorption   = cdom/nap absorption spectra, format = vector
%     wavelength   = wavelength values affiliated with absorption spectra,
%                    format = vector
%     cdom_model   = exponential, exponential_k, stretched_exponential or
%                    hyperbolic
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

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

%%

% Index matching input wavelength to specified minimum and maximum
% wavelength
ind=wavelength>=lambda_start & wavelength<=lambda_stop;

absorption=absorption(:);
wavelength=wavelength(:);

x=wavelength(ind);
y=absorption(ind);

ind=~isnan(y);
x=x(ind);
y=y(ind);

switch cdom_model
    case 'exponential'
        ft=fittype('cdom_model_exponential(x,a,s,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s'},'problem',{'lam0'}); % to make something a constant, specify it as 'problem',{'variable'}
        fopts=fitoptions(ft);
        fopts.StartPoint=[max(y) 0.01];
        [model,model_gof,model_fitstats]=fit(x,y,ft,fopts,'problem',{lam0}); % define problem variables, if any at the end as 'problem',{value1,value2}
    case 'exponential_k'
        ft=fittype('cdom_model_exponential_k(x,a,s,lam0,k)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','k'},'problem',{'lam0'});
        fopts=fitoptions(ft);
        fopts.StartPoint=[max(y) 0.01 mean(y(end-10:end))];
        fopts.Lower=[0 0 0];
        fopts.Upper=[inf 1 inf];
        [model,model_gof,model_fitstats]=fit(x,y,ft,fopts,'problem',{min(x)});
    case 'stretched_exponential'
        ft=fittype('cdom_model_stretched_exponential(x,a,s,beta,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','beta'},'problem',{'lam0'});
        fopts=fitoptions(ft);
        fopts.StartPoint=[max(y) 0.01 0.5];
        fopts.Lower=[0 0 0];
        fopts.Upper=[inf 1 1];
        [model,model_gof,model_fitstats]=fit(x,y,ft,fopts,'problem',{min(x)});
    case 'hyperbolic'
        ft=fittype('cdom_model_hyperbolic(x,a,s)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s'});
        fopts=fitoptions(ft);
        fopts.StartPoint=[max(y) 0.01];
        fopts.Lower=[0 0];
        fopts.Upper=[inf inf];
        [model,model_gof,model_fitstats]=fit(x,y,ft,fopts);
end


