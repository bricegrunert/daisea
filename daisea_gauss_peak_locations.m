function res=daisea_gauss_peak_locations(wavelength,phy)

% res=daisea_gauss_peak_locations(wavelength,phy)
%
% Find peak locations and number of peaks for gaussians.
%
% Inputs:
%  wavelength = wavelength
%  phy        = estimate of aph
%

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

% Created on 2022-06-10 by AC
% Steps 6 and 7 from DAISEA

%%
res.ff=[];
res.nphy=[];
res.number_gauss=0;
res.phi=[];
res.mu=[];
res.sigma=[];

% Ensure proper vector orientation
phy=phy(:);
wavelength=wavelength(:);

% Find wavelength interval from raw data
lam_int=unique(diff(wavelength));
if length(lam_int)>1
    error('Wavelength interval must be consistent.');
end

wavelength=wavelength(~isnan(phy));
phy=phy(~isnan(phy));

if length(phy)<5
    return;
end

ff=9/lam_int;
ff=round(ff,0);
if ff==2
    ff=3;
end

%nphy=sgolayfilt(phy,1,ff);
nphy=phy;

fn={'sigma1','phi1','mu1';'sigma2','phi2','mu2';'sigma3','phi3','mu3';'sigma4','phi4','mu4';'sigma5','phi5','mu5';'sigma6','phi6','mu6';'sigma7','phi7','mu7';'sigma8','phi8','mu8';'sigma9','phi9','mu9';'sigma10','phi10','mu10';'sigma11','phi11','mu11';'sigma12','phi12','mu12';'sigma13','phi13','mu13';'sigma14','phi14','mu14';'sigma15','phi15','mu15';'sigma16','phi16','mu16'};

drv_lam=[min(wavelength)+.5*lam_int:lam_int:max(wavelength)-.5*lam_int];
drv2_lam=[min(wavelength)+lam_int:lam_int:max(wavelength)-lam_int];

for jj=1:length(nphy)-1
    phy_drv(jj)=(nphy(jj)-nphy(jj+1))./lam_int;   %Tsai & Philpot 1998 Eq. 7
end

for jj=1:length(nphy)-2
    phy_drv2(jj)=(nphy(jj)-2*nphy(jj+1)+nphy(jj+2))/lam_int^2; %Tsai & Philpot 1998 Eq. 8
end

phy2=sgolayfilt(-phy_drv2,1,ff);
[phi,mu,sigma]=findpeaks(phy2,drv2_lam,'MinPeakHeight',0);

ind=find(sigma > 5);
phi=phi(ind);
mu=mu(ind);
sigma=sigma(ind);

[Y,I]=sort(phi,'descend');

[c,ia,ib]=intersect(mu,wavelength);

phy_res=nphy;

% Cap at 90% of absorption?
for ii=1:length(I)
    phi(I(ii))=phy_res(ib(I(ii))).*0.9;
    gg=gauss(wavelength,sigma(I(ii)),phi(I(ii)),mu(I(ii)));
    phy_res=phy_res-gg;
end

ind=find(phi > 0);

sigma=sigma(ind);
mu=mu(ind);
phi=phi(ind);

if length(phi) > 16
    [Y,I]=sort(phi,'descend');
    phi=phi(I(1:16));
    mu=mu(I(1:16));
    sigma=sigma(I(1:16));
end
res.mu1=mu;
res.phi1=phi;
res.sigma1=sigma;

if true
    % check gaussian fit results against pigment expectations
    % Nudge if >10, not cut
    % Minimum number of pigments - a&b
    pcheck.pigments={'Chl a','Chl a&c','Chl a','Chl b&c','Chl b','PPC','PSC','PE','Chl c','Chl a','Chl c','Chl b','Chl a'};
    pcheck.center=[382 409 436 457 468 491 525 552 585 620 639 658 676];
    pcheck.required=[true true true true true false false false false true false true true];
    %pcheck.required=[true true true false false false false false false true false false true];
    pcheck.halfwidth=[43 16 14 14 14 16 15 15 15 16 15 13 13];
    for ii=1:length(pcheck.center)
        if ~pcheck.required(ii); continue; end
        if pcheck.center(ii)<=min(wavelength); continue; end
        if pcheck.center(ii)>=max(wavelength); continue; end
        [tst,ind2]=min(abs(pcheck.center(ii)-mu));
        if tst>=10
            [tst2,ind3]=min(abs(wavelength-pcheck.center(ii)));
            mu=[mu pcheck.center(ii)];
            phi=[phi nphy(ind3)];
            sigma=[sigma pcheck.halfwidth(ii)*2];
        end
    end
    [junk,ind]=sort(mu);
    mu=mu(ind);
    phi=phi(ind);
    sigma=sigma(ind);
    phi(phi<0)=min(phi(phi>0));

    if length(phi) > 16
        [Y,I]=sort(phi,'descend');
        phi=phi(I(1:16));
        mu=mu(I(1:16));
        sigma=sigma(I(1:16));
    end
end

res.ff=ff;
res.nphy=nphy;
res.number_gauss=length(phi);
res.phi=phi;
res.mu=mu;
res.sigma=sigma;
