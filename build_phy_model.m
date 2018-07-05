function [model, model_gof, model_fitstats]=build_phy_model(absorption,wavelength,model,ngauss)

% build_phy_model
%
% [model, model_gof, model_fitstats]=build_phy_model(absorption,wavelength,model,ngauss)
%
% Function optimizes Gaussian curves previously identified for a given
% signal using a least squares fitting approach. For example, Gaussian curves identified on a phytoplankton
% absorption spectra using derivative analysis and findpeaks().
%
% Inputs:
%     absorption = input absorption spectra (or other signal) used to
%                  optimize input Gaussian components from a least squares
%                  approach, format = vector
%     wavelength = wavelength values affiliated with absorption spectra,
%                  format = vector
%     model      = model object with initial Gaussian estimates
%     ngauss     = number of Gaussian components used in model
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

fn3={'sigma1','phi1','mu1';'sigma2','phi2','mu2';'sigma3','phi3','mu3';'sigma4','phi4','mu4';'sigma5','phi5','mu5';'sigma6','phi6','mu6';'sigma7','phi7','mu7';'sigma8','phi8','mu8';'sigma9','phi9','mu9';'sigma10','phi10','mu10';'sigma11','phi11','mu11';'sigma12','phi12','mu12';'sigma13','phi13','mu13';'sigma14','phi14','mu14';'sigma15','phi15','mu15';'sigma16','phi16','mu16'};

if size(absorption,2)>1
    absorption=absorption';
end

if size(wavelength,2)>1
    wavelength=wavelength';
end

x=wavelength;
y=absorption;

ind=find(isnan(y)==0);

x=x(ind);
y=y(ind);

for ii=1:ngauss
    for jj=1:3
        nparam(ii,jj)=model.model_fit.(fn3{ii,jj});
    end
end
    
    
    wm=['phy_model_' num2str(ngauss) 'gaussian'];
    
    
    switch wm
        
        case 'phy_model_1gaussian'
            ft=fittype('gauss1(x,sigma1,phi1,mu1)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1'});% to make something a constant, specify it as 'problem',{'variable'}
            
        case 'phy_model_2gaussian'
            ft=fittype('gauss2(x,sigma1,phi1,mu1,sigma2,phi2,mu2)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2'});
            
        case 'phy_model_3gaussian'
            ft=fittype('gauss3(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3'});
            
        case 'phy_model_4gaussian'
            ft=fittype('gauss4(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4'});
            
        case 'phy_model_5gaussian'
            ft=fittype('gauss5(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5'});
            
        case 'phy_model_6gaussian'
            ft=fittype('gauss6(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6'});
            
        case 'phy_model_7gaussian'
            ft=fittype('gauss7(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7'});
            
        case 'phy_model_8gaussian'
            ft=fittype('gauss8(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8'});
            
        case 'phy_model_9gaussian'
            ft=fittype('gauss9(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9'});
            
        case 'phy_model_10gaussian'
            ft=fittype('gauss10(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10'});
            
        case 'phy_model_11gaussian'
            ft=fittype('gauss11(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11'});
            
        case 'phy_model_12gaussian'
            ft=fittype('gauss12(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12'});
            
        case 'phy_model_13gaussian'
            ft=fittype('gauss13(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12','sigma13','phi13','mu13'});
            
        case 'phy_model_14gaussian'
            ft=fittype('gauss14(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13,sigma14,phi14,mu14)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12','sigma13','phi13','mu13','sigma14','phi14','mu14'});
            
        case 'phy_model_15gaussian'
            ft=fittype('gauss15(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13,sigma14,phi14,mu14,sigma15,phi15,mu15)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12','sigma13','phi13','mu13','sigma14','phi14','mu14','sigma15','phi15','mu15'});

            case 'phy_model_16gaussian'
            ft=fittype('gauss16(x,sigma1,phi1,mu1,sigma2,phi2,mu2,sigma3,phi3,mu3,sigma4,phi4,mu4,sigma5,phi5,mu5,sigma6,phi6,mu6,sigma7,phi7,mu7,sigma8,phi8,mu8,sigma9,phi9,mu9,sigma10,phi10,mu10,sigma11,phi11,mu11,sigma12,phi12,mu12,sigma13,phi13,mu13,sigma14,phi14,mu14,sigma15,phi15,mu15,sigma16,phi16,mu16)','independent',{'x'},'dependent',{'y'},'coefficients',{'sigma1','phi1','mu1','sigma2','phi2','mu2','sigma3','phi3','mu3','sigma4','phi4','mu4','sigma5','phi5','mu5','sigma6','phi6','mu6','sigma7','phi7','mu7','sigma8','phi8','mu8','sigma9','phi9','mu9','sigma10','phi10','mu10','sigma11','phi11','mu11','sigma12','phi12','mu12','sigma13','phi13','mu13','sigma14','phi14','mu14','sigma15','phi15','mu15','sigma16','phi16','mu16'});

        otherwise
            error(['No match for model - time to make a new one! ' wm]);
    end
        
        
        fn={'mu1','mu2','mu3','mu4','mu5','mu6','mu7','mu8','mu9','mu10','mu11','mu12','mu13','mu14','mu15','mu16'};
        fn2={'sigma1','sigma2','sigma3','sigma4','sigma5','sigma6','sigma7','sigma8','sigma9','sigma10','sigma11','sigma12','sigma13','sigma14','sigma15','sigma16'};
        fn4={'phi1','phi2','phi3','phi4','phi5','phi6','phi7','phi8','phi9','phi10','phi11','phi12','phi13','phi14','phi15','phi16'};

        fopts=fitoptions(ft);
        for cc=1:ngauss
            fopts.StartPoint=[fopts.StartPoint model.model_fit.(fn2{cc}) model.model_fit.(fn4{cc}) model.model_fit.(fn{cc})];
        end        
        
        for cc=1:ngauss
            fopts.Upper=[fopts.Upper model.model_fit.(fn2{cc})*2 model.model_fit.(fn4{cc})*1.5 model.model_fit.(fn{cc})+5];
            fopts.Lower=[fopts.Lower model.model_fit.(fn2{cc})*.5 model.model_fit.(fn4{cc})*.25 model.model_fit.(fn{cc})-5];
        end

    
    [model,model_gof,model_fitstats]=fit(x,y,ft,fopts);
    
