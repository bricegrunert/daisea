function [model, model_gof, model_fitstats]=build_daisea_model(absorption,wavelength,model,ngauss,Sdg_est,lam0,alam0,lambda_start,lambda_stop)

% build_dg_model
%
% [model, model_gof, model_fitstats]=build_daisea_model(absorption,wavelength,lambda_start,lambda_stop,lam0,model,ngauss,peak_lim,Sdg_est,alam0)
%
% Function identifies a suitable model to fit total non-water absorption
% spectra using a single exponential and pre-determined number of Gaussian 
% curves. Steps to initially estimate input values derived from steps in
% daisea.m
%
% Inputs:
%     absorption   = input absorption spectra (or other signal) used to
%                    optimize input Gaussian components from a least
%                    squares approach, format = vector
%     wavelength   = wavelength values affiliated with absorption spectra,
%                    format = vector
%     model        = model object with initial Gaussian estimates
%     ngauss       = number of Gaussian components used in model
%     Sdg_est      = input estimate of NAP/CDOM absorption spectral slope
%     lam0         = wavelength used to initialize fitting of NAP/CDOM
%                    absorption
%     alam0        = input estimate of NAP/CDOM absorption at lam0
%     lambda_start = minimum wavelength considered (e.g. 350 nm),
%                    must be within spectral range of wavelength input
%     lambda_stop  = maximum wavelength considered (e.g. 700 nm),
%                    must be within spectral range of wavelength input
%
% Returns:
%     model          = output from Gaussian decomposition of
%                      absorption
%     model_gof      = goodness of fit for model
%     model_fitstats = fitting stats for model
% 
% copyright (c) 2018 Brice K. Grunert
% email: bricegrunert@gmail.com

%Last modified on 28 June 2018 by BG


%%

if nargin==7
    lambda_start=min(wavelength);
    lambda_stop=max(wavelength);
end

fn3={'sigma1','phi1','mu1';'sigma2','phi2','mu2';'sigma3','phi3','mu3';'sigma4','phi4','mu4';'sigma5','phi5','mu5';'sigma6','phi6','mu6';'sigma7','phi7','mu7';'sigma8','phi8','mu8';'sigma9','phi9','mu9';'sigma10','phi10','mu10';'sigma11','phi11','mu11';'sigma12','phi12','mu12';'sigma13','phi13','mu13';'sigma14','phi14','mu14';'sigma15','phi15','mu15';'sigma16','phi16','mu16'};

for jj=1:length(wavelength)
    if wavelength(jj)==lambda_start
        ind=jj;
    end
    if wavelength(jj)==lambda_stop
        ind2=jj;
    end
end

lind=find(wavelength==lam0);

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

for ii=1:ngauss
    for jj=1:3
        nparam(ii,jj)=model.model_fit.(fn3{ii,jj});
    end
end

wm=['cdom_model_' num2str(ngauss) 'gaussian'];

switch wm
    
    case 'cdom_model_0gaussian'
        ft=fittype('cdom_model_0gaussian_noK(x,a,s,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s'},'problem',{'lam0'});% to make something a constant, specify it as 'problem',{'variable'}
        
    case 'cdom_model_1gaussian'
        ft=fittype('cdom_model_1gaussian_noK(x,a,s,sigma1,phi1,mu1,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1'},'problem',{'lam0'});% to make something a constant, specify it as 'problem',{'variable'}
        
    case 'cdom_model_2gaussian'
        ft=fittype('cdom_model_2gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2'},'problem',{'lam0'});
        
    case 'cdom_model_3gaussian'
        ft=fittype('cdom_model_3gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3'},'problem',{'lam0'});
        
    case 'cdom_model_4gaussian'
        ft=fittype('cdom_model_4gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4'},'problem',{'lam0'});
        
    case 'cdom_model_5gaussian'
        ft=fittype('cdom_model_5gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5'},'problem',{'lam0'});
        
    case 'cdom_model_6gaussian'
        ft=fittype('cdom_model_6gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6'},'problem',{'lam0'});
        
    case 'cdom_model_7gaussian'
        ft=fittype('cdom_model_7gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7'},'problem',{'lam0'});
        
    case 'cdom_model_8gaussian'
        ft=fittype('cdom_model_8gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8'},'problem',{'lam0'});
        
    case 'cdom_model_9gaussian'
        ft=fittype('cdom_model_9gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9'},'problem',{'lam0'});
        
    case 'cdom_model_10gaussian'
        ft=fittype('cdom_model_10gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10'},'problem',{'lam0'});
        
    case 'cdom_model_11gaussian'
        ft=fittype('cdom_model_11gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11'},'problem',{'lam0'});
        
    case 'cdom_model_12gaussian'
        ft=fittype('cdom_model_12gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12'},'problem',{'lam0'});
        
    case 'cdom_model_13gaussian'
        ft=fittype('cdom_model_13gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12','sigma13','phi13','mu13'},'problem',{'lam0'});
        
    case 'cdom_model_14gaussian'
        ft=fittype('cdom_model_14gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13,sigma14,phi14,mu14,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12','sigma13','phi13','mu13','sigma14','phi14','mu14'},'problem',{'lam0'});
        
    case 'cdom_model_15gaussian'
        ft=fittype('cdom_model_15gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13,sigma14,phi14,mu14,sigma15,phi15,mu15,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12','sigma13','phi13','mu13','sigma14','phi14','mu14','sigma15','phi15','mu15'},'problem',{'lam0'});
        
    case 'cdom_model_16gaussian'
        ft=fittype('cdom_model_16gaussian_noK(x,a,s,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13,sigma14,phi14,mu14,sigma15,phi15,mu15,sigma16,phi16,mu16,lam0)','independent',{'x'},'dependent',{'y'},'coefficients',{'a','s','sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12','sigma13','phi13','mu13','sigma14','phi14','mu14','sigma15','phi15','mu15','sigma16','phi16','mu16'},'problem',{'lam0'});
        
    otherwise
        error(['No match for model - time to make a new one! ' wm]);
end

if ngauss > 0
    lr=size(nparam,1);
    lc=size(nparam,2);
    
    nparam=reshape(nparam',lr*lc,1);
    
    fopts=fitoptions(ft);
    fopts.StartPoint=[alam0 Sdg_est nparam'];
    
    fn={'mu1','mu2','mu3','mu4','mu5','mu6','mu7','mu8','mu9','mu10','mu11','mu12','mu13','mu14','mu15','mu16'};
    fn2={'sigma1','sigma2','sigma3','sigma4','sigma5','sigma6','sigma7','sigma8','sigma9','sigma10','sigma11','sigma12','sigma13','sigma14','sigma15','sigma16'};
    fn4={'phi1','phi2','phi3','phi4','phi5','phi6','phi7','phi8','phi9','phi10','phi11','phi12','phi13','phi14','phi15','phi16'};
    
    fopts.Upper=[absorption(lind) Sdg_est+0.003];
    fopts.Lower=[0 Sdg_est-0.002];
    
    for cc=1:ngauss
        fopts.Upper=[fopts.Upper model.model_fit.(fn2{cc})*3 model.model_fit.(fn4{cc})*3 model.model_fit.(fn{cc})];
        fopts.Lower=[fopts.Lower model.model_fit.(fn2{cc}) model.model_fit.(fn4{cc})*.25 model.model_fit.(fn{cc})];
    end
    
else
    fopts=fitoptions(ft);
    fopts.StartPoint=[alam0 Sdg_est];
    fopts.Upper=[absorption(lind) Sdg_est+0.003];
    fopts.Lower=[0 Sdg_est-0.002];
end

[model, model_gof, model_fitstats]=fit(x,y,ft,fopts,'problem',[lam0]);
