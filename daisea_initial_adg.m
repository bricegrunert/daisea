function [output,qcflag,ind]=daisea_initial_adg(output,wlind,cdom_model)

% [output,qcflag]=daisea_initial_adg(output,wlind)
%
% Steps 1 and 2 - Estimate initial adg. Calcualte second derivative, use
% selected points to retrieve slope. 
%
% Inputs:
%  output     = output structrue from daisea
%  wlind      = index structure for wavelength positions
%  cdom_model = exponential, exponential_k, stretched_exponential,
%               hyperbolic
%

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

% Created on 2022-09-22 by AC

%% Step 1
qcflag=false;

% calculate 1st and 2nd derivative
% finite approximation derivative
output.drv=(output.anw(1:end-1)-output.anw(2:end))./output.wavelength_interval; %Tsai & Philpot 1998 Eq. 7
output.drv2=(output.anw(1:end-2)-2*output.anw(2:end-1)+output.anw(3:end))./output.wavelength_interval^2; %Tsai & Philpot 1998 Eq. 8

output.drv=output.drv';
output.drv2=output.drv2';

tst=std(output.drv2);
ind=find(output.drv2<=tst & output.drv2>=-tst);

% exclude red chl peak, 676 +/- 15  nm
[junk,lb]=min(abs(output.wavelength-(676-15)));
[junk,ub]=min(abs(output.wavelength-(676+15)));
ind=ind(ind<lb | ind>ub);
% exclude blue chl peak, 457 +/- 15  nm
[junk,lb]=min(abs(output.wavelength-(457-15)));
[junk,ub]=min(abs(output.wavelength-(457+15)));
ind=ind(ind<lb | ind>ub);

% if there are less than 5 points where the 2nd deriv = 0, do
% not attempt to decompose the spectra
if length(ind)<5
    qcflag=true;
    return;
end

% if more than 4 points where 2nd deriv=0 are found, proceed
% calculate spectral slope on total absorption where 2nd deriv=0
output.n_points=length(ind);
lam_start=output.wavelength(ind(1));
lam_stop=output.wavelength(ind(end));

%% Step 2
[Sg,Sg_g,Sg_o]=daisea_fit_cdom_nap(output.anw(ind),output.wavelength(ind),cdom_model,lam_start,lam_stop,440);

%Initial estimates
output.adg_model1=daisea_apply_cdom_model(cdom_model,output.wavelength,'fit',Sg);
output.adg_fit1=Sg;
output.adg_model1_adjusted=daisea_apply_cdom_model(cdom_model,output.wavelength,'fit',Sg,'percent_dg',output.percent_dg);

