function [model, model_gof, model_fitstats]=daisea_gaussian_fit(absorption,wavelength,gaussian_fit,ngauss)

% [model, model_gof, model_fitstats]=daisea_gaussian_fit(absorption,wavelength,gaussian_fit,ngauss)
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
%     gaussian_fit = model object with initial Gaussian estimates
%     ngauss     = number of Gaussian components used in model
%
% Returns:
%     model          = output from Gaussian decomposition of
%                      absorption
%     model_gof      = goodness of fit for model
%     model_fitstats = fitting stats for model
%

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

% Created on 2022-12-15 by AC

%%
absorption=absorption(:);
wavelength=wavelength(:);

x=wavelength(~isnan(absorption));
y=absorption(~isnan(absorption));

wm=['gauss' num2str(ngauss)];
nm=[wm '(x'];

fn={'sigma','phi','mu'};
sp=NaN(1,ngauss*length(fn));
c=cell(size(sp));
cnt=0;
for ii=1:ngauss
    for jj=1:length(fn)
        cnt=cnt+1;
        sp(cnt)=gaussian_fit.([fn{jj} num2str(ii)]);
        c{cnt}=[fn{jj} num2str(ii)];
        nm=[nm ',' [fn{jj} num2str(ii)]];
    end
end
nm=[nm ')'];

ft=fittype(nm,'independent',{'x'},'dependent',{'y'},'coefficients',c);
fopts=fitoptions(ft);
fopts.StartPoint=sp;
for ii=1:ngauss
    for jj=1:length(fn)
        switch fn{jj}
            case 'sigma'
                fopts.Upper=[fopts.Upper gaussian_fit.([fn{jj} num2str(ii)])*2];
                fopts.Lower=[fopts.Lower gaussian_fit.([fn{jj} num2str(ii)])*0.5];
                %fopts.Upper=[fopts.Upper max(wavelength)-min(wavelength)];
                %fopts.Lower=[fopts.Lower 1];
            case 'phi'
                fopts.Upper=[fopts.Upper gaussian_fit.([fn{jj} num2str(ii)])*1.5];
                fopts.Lower=[fopts.Lower gaussian_fit.([fn{jj} num2str(ii)])*0.25];
                %fopts.Upper=[fopts.Upper max(absorption)];
                %fopts.Lower=[fopts.Lower 0];
            case 'mu'
                fopts.Upper=[fopts.Upper gaussian_fit.([fn{jj} num2str(ii)])+10];
                fopts.Lower=[fopts.Lower gaussian_fit.([fn{jj} num2str(ii)])-10];
        end
    end
end

[model,model_gof,model_fitstats]=fit(x,y,ft,fopts);

