function output=daisea(absorption,wavelength,settings)

% daisea
%
% output=daisea(absorption,wavelength)
%
% Function uses measured total absorption data and decomposes into
% phytoplankton and NAP/CDOM components using pre-filters and initial
% estimates to constrain a Gaussian decomposition of total non-water 
% absorption using least squares fitting. 
%
% output=daisea(absorption,wavelength,settings) allows the user to provide
% settings for data processing. See description of "Inputs" below. 
%
% Approach is intended for hyperspectral data from 350-700 nm, with a
% maximum bandwidth of 5 nm. Algorithm can handle data at a finer
% resolution. Algorithm specifically uses 350, 440, 555, 680 and 690 nm
% wavelengths but will find values near these by selecting the nearest
% wavelength. If input wavelengths are significantly different from these
% values, algorithm performance will likely be poor.
%
% Users should verify ability to perform beyond this spectral window
% and resolution before use. 
%
% Steps below refer to manuscript:
% Deriving inherent optical properties from decomposition of hyperspectral
% non-water absorption
% B.K. Grunert, C.B. Mouw, A.B. Ciochetto
%
% We highly recommend users refer to this paper for a detailed description
% of each step with an accompanying flow chart as well as providing a
% thorough description of model output.
%
% Inputs:
%     absorption = total non-water absorption spectra (format = vector)
%     wavelength = wavelength values affiliated with absorption spectra
%                  (format = vector)
%     settings   = user defined settings (format = structure), optional, if
%                  not provided by user, defaults will be applied as
%                  descirbed below. 
%
% Fields for settings structure:
%    Note that all fields are optional. If a field is excluded from the
%    structure, it will be set to the default. If a field is set to the
%    empty set or NaN, it will also be set to the default. The user can
%    change the defaults at the start of the code below. 
% 
%    notify             = true or false, print warnings to screen, default
%                         is true
%    minimum_wavelength = minimum wavelength considered (e.g. 300 nm),
%                         default is min(wavelength)
%    maximum_wavelength = maximum wavelength considered (e.g. 700 nm),
%                         default is max(wavelength)
%    negative_data      = DAISEA requires all data to be positive. Negative
%                         data can either be 'offset' (add lowest value to
%                         entire dataset to make it positive), 'zero' (set 
%                         to 0), or 'remove' (set to NaN and ignored),
%                         default is 'offset'
%    wavelength_interval = DAISEA requires all data to have a regular 
%                         wavelength interval. If they do not (e.g. ACS 
%                         data), then, by default the program will 
%                         interpolate the data to a wavelength interval of
%                         floor(median(diff(wavelength))). Otherwise, the 
%                         user may enter a value here to override the 
%                         default. 
% 
%                         If your data have a regular wavelength 
%                         interval and you enter a value here that is 
%                         different from that interval, then the data will
%                         be interpolated to the interval specified here. 
%
%                         If your data have a regular wavelength
%                         interval but get removed (e.g. because they are
%                         negative or NaN), then the truncated data set
%                         will be interpolated back to the original
%                         wavelength interval. Again, if you enter a value
%                         here that is different from the wavelength
%                         interval of your data, then the data will be
%                         interpolated to the interval specified here. 
%
% Returns:
%     output = structure containing information from derivative analysis 
%     and iterative spectral fitting or from full DAISEA algorithm. If a
%     field is unspecified, it is left blank in the structure. Fields are
%     as follows:
%
%           at_nw                = input total non-water absorption
%           lam                  = input wavelengths for at_nw
%           drv                  = first derivative of total non-water
%                                  absorption
%           drv2                 = second derivative of total non-water
%                                  absorption
%           n_points             = number of data points where drv2 = 0
%           initial_adg_model    = adg estimated using exponential fit of
%                                  at_nw at spectral locations of n_points
%           Sgd_initial_est      = Sgd estimate for initial_adg_model
%           adg_model2           = second estimated adg after spectral 
%                                  evaluation 
%           model2_agd_lam0      = estimated agd at lam0 for adg_model2
%           model2_lam0          = lam0 used for adg_model2
%           Sgd_est2             = Sgd estimated for adg_model2
%                                 relaxation
%           final_model          = output from Gaussian decomposition of
%                                  at_nw
%           final_model_gof      = goodness of fit for final_model
%           final_model_fitstats = fitting stats for final_model
%           agd_final_estimate   = length of multiple turnover relaxation
%           Sgd_final_est        = Sgd estimate from agd_final_estimate
%           aphy_final_estimate  = estimate of aphy from at_nw -
%                                  agd_final_estimate
%           dg_dominant          = 'Yes' or 'No' estimating whether percent
%                                  contribution of aphy(440) falls below 
%                                  10% 
%           final_lam0           = lam0 used for agd_final_estimate
%           agd_final_lam0       = agd_final_estimate value at final_lam0
%           number_gauss         = number of Gaussian components used to
%                                  parameterize aphy in final_model 
% 
% copyright (c) 2018 Brice K. Grunert
% email: bricegrunert@gmail.com

%Last modified on 20 June 2018 by BG


%% Step 0
% Re-assign input absorption and wavelength to appear in output structure in
% the proper format (single column for each parameter), remove NaN's and
% negative values, assign start and stop lambda if needed

% Create settings structure and check on user input
if nargin==2 % defaults for settings structure
    settings.notify              = true;
    settings.minimum_wavelength  = min(wavelength);
    settings.maximum_wavelength  = max(wavelength);
    settings.negative_data       = 'offset'; % offset, zero or remove
    settings.wavelength_interval = []; % allows user to specify wavelength interval for interpoltion of unevenly spaced data
else
    if ~isfield(settings,'minimum_wavelength') || isempty(settings.minimum_wavelength) || isnan(settings.minimum_wavelength) || ischar(settings.minimum_wavelength)
        settings.minimum_wavelength  = min(wavelength);
    end
    if ~isfield(settings,'maximum_wavelength') || isempty(settings.maximum_wavelength) || isnan(settings.maximum_wavelength) || ischar(settings.maximum_wavelength)
        settings.maximum_wavelength  = max(wavelength);
    end
    if ~isfield(settings,'negative_data') || isnumeric(settings.negative_data)
        settings.negative_data='offset';
    end
    if ~isfield(settings,'wavelength_interval') || ischar(settings.wavelength_interval)
        settings.wavelength_interval=[];
    end
    if ~isfield(settings,'notify') || ischar(settings.wavelength_interval)
        settings.notify=true;
    end
end

% Initialize structure
output.at_nw=[];
output.lam=[];
output.drv=[];
output.drv2=[];
output.n_points=[];
output.initial_adg_model=[];
output.Sgd_initial_est=[];
output.adg_model2=[];
output.model2_agd_lam0=[];
output.model2_lam0=[];
output.Sgd_est2=[];
output.final_model=[];
output.final_model_gof=[];
output.final_model_fitstats=[];
output.agd_final_estimate=[];
output.Sgd_final_est=[];
output.aphy_final_estimate=[];
output.dg_dominant=[];
output.final_lam0=[];
output.agd_final_lam0=[];
output.number_gauss=[];

% Ensure proper vector orientation
absorption=absorption(:);
wavelength=wavelength(:);

% Find wavelength interval from raw data
lam_int=unique(diff(wavelength));

% Get rid of NaNs
wavelength=wavelength(~isnan(absorption));
absorption=absorption(~isnan(absorption));
if isempty(absorption)
    error('ERROR DAISEA: Data are all NaN.');
end

% If minimum is less than 0, add offset or set to zero
if min(absorption)<0
    switch settings.negative_data
        case 'offset'
            absorption=absorption-min(absorption);
        case 'zero'
            absorption(absorption<0)=0;
        case 'remove'
            wavelength=wavelength(absorption>0);
            absorption=absorption(absorption>0);
        otherwise
            error(['ERROR DAISEA: ' settings.negative_data ' is not an option for settings.negative_data. Enter ''offset'', ''zero'', or ''remove''.']);
    end
end

% If wavelength interval is not consistent, fix
lam_int2=unique(diff(wavelength));
fixme=false;
if (length(lam_int)==1 && length(lam_int2)>1)
    fixme=true; % fix but use original wavelength interval
elseif length(lam_int)>1 && length(lam_int2)>1
    fixme=true; % fix and determine wavelength interval
    lam_int=floor(median(diff(wavelength)));
end
% get wavelength interval
if any((~isnan(settings.wavelength_interval) | ~isempty(settings.wavelength_interval)) & settings.wavelength_interval~=lam_int)
    fixme=true;
    lam_int=settings.wavelength_interval;
end
if lam_int>5
    fprintf('\nWARNING DAISEA: Wavelength interval is %0.2f nm. Data will be interpolated to the maximum allowable interval of 5 nm.',lam_int);
    fixme=true;
    lam_int=5;
end
% Put onto new grid
if fixme
    if settings.notify
        fprintf('\nDAISEA: Data are being interpolated to a %0.2f nm wavelength interval',lam_int);
    end
    wavelength2=ceil(min(wavelength)):lam_int:floor(max(wavelength));
    absorption2=interp1(wavelength,absorption,wavelength2,'linear',NaN);
    output.lam=wavelength2(:);
    output.at_nw=absorption2(:);
else
    output.lam=wavelength;
    output.at_nw=absorption;
end
clear lam_int2 fixme wavelength2 absorption2

%Find data falling within the desired fitting window (e.g. 350-700 nm)
[start,start_ind]=min(abs(output.lam-settings.minimum_wavelength));
[stop,stop_ind]=min(abs(output.lam-settings.maximum_wavelength));

start=output.lam(start_ind);
stop=output.lam(stop_ind);

output.lam=output.lam(start_ind:stop_ind);
output.at_nw=output.at_nw(start_ind:stop_ind);

%In situ data - calculate difference in total absorption spectral slope for
%all wavelengths and wavelengths outside the chl-a peak region (410-480 nm)
[junk,ind555]=min(abs(output.lam-555));
if junk>5 && settings.notify
    fprintf('\nWARNING DAISEA: Closest wavlenth to 555 nm is %0.2f nm away.',junk);
end
[junk,ind680]=min(abs(output.lam-680));
if junk>5 && settings.notify
    fprintf('\nWARNING DAISEA: Closest wavlenth to 680 nm is %0.2f nm away.',junk);
end
[junk,ind690]=min(abs(output.lam-690));
if junk>5 && settings.notify
    fprintf('\nWARNING DAISEA: Closest wavlenth to 690 nm is %0.2f nm away.',junk);
end
[junk,ind350]=min(abs(output.lam-350));
if junk>5 && settings.notify
    fprintf('\nWARNING DAISEA: Closest wavlenth to 350 nm is %0.2f nm away.',junk);
end
[junk,ind440]=min(abs(output.lam-440));
if junk>5 && settings.notify
    fprintf('\nWARNING DAISEA: Closest wavlenth to 440 nm is %0.2f nm away.',junk);
end
[junk,ind400]=min(abs(output.lam-400));
if junk>5 && settings.notify
    fprintf('\nWARNING DAISEA: Closest wavlenth to 400 nm is %0.2f nm away.',junk);
end
clear junk
if ind400==ind350 && settings.notify
    fprintf('\nWARNING DAISEA: Closest wavelength to both 350 nm and 400 nm is %0.2f nm',output.lam(ind400));
end
    
ratio=output.at_nw(ind555)/output.at_nw(ind680);

%if the slope difference is greater than 0.0003, treat as a
%predominantly NAP/CDOM sample
%if ratio > 3.07
if ratio > 2.528
    output.dg_dominant='Yes';
    
    %Initial and final estimates are the same, calcualted from total
    %water absorption (more accurate than finding where 2nd deriv=0 and
    %calculating for that)
    
    [Sg, Sg_g, Sg_o]=spectral_slope_nogauss_noK(output.at_nw,output.lam,start,stop,start);
    output.Sgd_initial_est=Sg.s;
    output.initial_adg_model=cdom_model_lam0_noK(output.lam,output.at_nw(ind350),output.Sgd_initial_est,start);
    output.Sgd_final_est=Sg.s;
    
    %are there wavelengths where estimated agd is greater than observed
    %at? find max overestimate and use this as the offset, with that
    %wavelength as the new lam0
    tst=output.initial_adg_model-output.at_nw;
    nind=find(tst > 0);
    lam_ind=find(tst==max(tst));
    offset=tst(lam_ind);
    
    %fit an exponential curve based off the new lam0 and offset value
    if output.at_nw(lam_ind)-offset > 0
        output.agd_final_estimate=cdom_model_lam0_noK(output.lam,output.at_nw(lam_ind)-offset,output.Sgd_final_est,output.lam(lam_ind));
        output.final_lam0=output.lam(lam_ind);
        output.agd_final_lam0=output.at_nw(lam_ind)-offset;
    else
        output.agd_final_estimate=output.initial_adg_model;
        output.final_lam0=start;
        output.agd_final_lam0=output.at_nw(ind350);
    end
    
    ff=9/lam_int;
    ff=round(ff,0);
    if ff==2
        ff=3;
    end
    
    output.aphy_final_estimate=sgolayfilt(output.at_nw-output.agd_final_estimate,1,ff);
    
    %Create empty variables for those used in decomposition below but
    %not here
    output.drv=[];
    output.drv2=[];
    output.n_points=[];
    output.model2_agd_lam0=[];
    output.model2_lam0=[];
    output.adg_model2=[];
    output.Sgd_est2=[];
    output.final_model=[];
    output.final_model_gof=[];
    output.final_model_fitstats=[];
    output.number_gauss=[];
    
    
else
    %Samples with a significant phytoplankton contribution
    output.dg_dominant='No';
    
    % empirical relationship from training dataset - Grunert et al.
    % submitted June 2018
    if output.at_nw(ind555)/output.at_nw(ind680) > 0.685
        phyperc=1.038*exp(-.9257*(output.at_nw(ind555)/output.at_nw(ind680)));
    else
        phyperc=2.088*exp(-1.946*(output.at_nw(ind555)/output.at_nw(ind680)));
    end
    dgperc=1-phyperc;
    
%% Step 1

    %calculate 1st and 2nd derivative
    %finite approximation derivative
    for jj=1:length(output.at_nw)-1
        output.drv(jj)=(output.at_nw(jj)-output.at_nw(jj+1))./lam_int;   %Tsai & Philpot 1998 Eq. 7
    end
    for jj=1:length(output.at_nw)-2
        output.drv2(jj)=(output.at_nw(jj)-2*output.at_nw(jj+1)+output.at_nw(jj+2))/lam_int^2; %Tsai & Philpot 1998 Eq. 8
    end
    
    mag=-floor(log10(nanmean(output.at_nw)))+2;
    if mag>=0
        ind=find(round(output.drv2,mag)==0);
        while length(ind) < 5 & mag >= 0
            mag=mag-1;
            ind=find(round(output.drv2,mag)==0);
        end
    else
        mag=0;
        ind=find(round(output.drv2,mag)==0);
    end
    
    %if more than 4 points where 2nd deriv=0 are found, proceed
    if length(ind) >= 5
        
        %calculate spectral slope on total absorption where 2nd deriv=0
        output.n_points=length(ind);
        lam_start=output.lam(ind(1));
        lam_stop=output.lam(ind(end));
        
%% Step 2

        [Sg, Sg_g, Sg_o]=spectral_slope_nogauss_noK(output.at_nw(ind),output.lam(ind),lam_start,lam_stop,lam_start);
        %Initial slope estimate
        output.Sgd_initial_est=Sg.s;
        
        %create vector to alter input slope
        itVec=[0:0.0001:0.011];
        itVec2=[-0.004:.0001:-0.0001];
        itVec(length(itVec)+1:length(itVec)+length(itVec2))=sort(itVec2,'descend');
        
        output.model2_agd_lam0=NaN;
        output.model2_lam0=NaN;
        output.Sgd_est2=NaN;
        output.initial_adg_model=cdom_model_lam0_noK(output.lam,output.at_nw(ind440)*dgperc,output.Sgd_initial_est,output.lam(ind440));
        
        if min(output.at_nw(end)) < 0
            tst=output.initial_adg_model(ind350:ind690)-output.at_nw(ind350:ind690);
        else
            tst=output.initial_adg_model-output.at_nw;
        end
        
        %are there points where estimated agd is above measured at_nw?
        nind=find(tst > 0);
        
        %if the original agd estimation is below at_nw at all points,
        %stick with it
        
%% Step 3

        if length(nind)==0
            output.Sgd_est2=output.Sgd_initial_est;
            output.model2_lam0=output.lam(ind440);
            output.model2_agd_lam0=output.at_nw(ind440)*dgperc;
            output.adg_model2=output.initial_adg_model;
        else
            
            off_ind=find(tst==max(tst));
            offset=tst(off_ind);

            %if original estimation is above, iterate through different
            %slope values until agd is less than at_nw at all
            %wavelengths
            cnt=0;
            itWorks=zeros(length(itVec),1);
            
            while cnt <= 150 && sum(itWorks)==0
                
                cnt=cnt+1;
                
                S_in=output.Sgd_initial_est+itVec(cnt);
                tst2(:,cnt)=cdom_model_lam0_noK(output.lam,output.at_nw(off_ind),S_in,output.lam(off_ind));
                est2(:,cnt)=tst2(1:end-1,cnt)-output.at_nw(1:end-1);
                if max(est2(:,cnt)) > 0
                    output.Sgd_est2=NaN;
                    output.model2_lam0=NaN;
                    output.model2_agd_lam0=NaN;
                    itWorks(cnt)=0;
                else
                    output.Sgd_est2=output.Sgd_initial_est+itVec(cnt);
                    output.model2_lam0=output.lam(ind440);
                    output.model2_agd_lam0=output.at_nw(ind440)*dgperc;
                    itWorks(cnt)=1;
                end
            end
            
            output.adg_model2=cdom_model_lam0_noK(output.lam,output.model2_agd_lam0,output.Sgd_est2,output.model2_lam0);
                        
        end

        %if we found a spectra that satisfies the above requirements,
        %proceed with the decomposition
        if isnan(output.Sgd_est2)==1
            
            output.Sgd_est2=output.Sgd_initial_est;
            output.model2_lam0=output.lam(ind440);
            output.model2_agd_lam0=output.at_nw(ind440)*dgperc;
            output.adg_model2=cdom_model_lam0_noK(output.lam,output.model2_agd_lam0,output.Sgd_est2,output.model2_lam0);
            
        end

%% Step 4

        phy=output.at_nw-output.adg_model2;
                
%% Step 5

        if phy(ind350)/phy(ind440) > 1.5
            
            S=spectral_slope_nogauss_noK(phy,output.lam,start,output.lam(ind400),start);
            yy=cdom_model_lam0_noK(output.lam,S.a,S.s,start);    
            nyy=output.adg_model2+yy;
            
            S=spectral_slope_nogauss_noK(nyy,output.lam,start,stop,start);            
            alam0=output.at_nw(ind440)*dgperc;
            lam_ind=ind440;
            agd_tst=cdom_model_lam0_noK(output.lam,alam0,S.s,output.lam(lam_ind));
            
            clear tst
            tst=agd_tst(ind350:ind690)-output.at_nw(ind350:ind690);
            
            if max(tst) > 0
                off_ind=find(tst==max(tst));
                offset=tst(off_ind);
                
                alam0=agd_tst(off_ind)-offset;
                lam_ind=off_ind;
                agd_tst=cdom_model_lam0_noK(output.lam,alam0,S.s,output.lam(off_ind));
            end
            
            phy=output.at_nw-agd_tst;
            output.Sgd_est2=S.s;
            output.adg_model2=agd_tst;
            output.model2_lam0=output.lam(lam_ind);
            output.model2_agd_lam0=alam0;
            
            cnt=0;
            
            while phy(ind350)/phy(ind440) > 1.5 & cnt < 20
                
                cnt=cnt+1;
                
                S_in=S.s+cnt*0.0001;
                
                agd_tst=cdom_model_lam0_noK(output.lam,output.model2_agd_lam0,S_in,output.model2_lam0);
                
                clear tst
                tst=agd_tst(ind350:ind690)-output.at_nw(ind350:ind690);
                
                if max(tst) > 0
                    lam_ind=find(tst==max(tst));
                    offset=tst(lam_ind);
                    alam0=output.at_nw(lam_ind)-offset;
                    agd_tst=cdom_model_lam0_noK(output.lam,alam0,S_in,output.lam(lam_ind));
                    output.model2_agd_lam0=alam0;
                    output.model2_lam0=output.lam(lam_ind);
                end
                
                phy=output.at_nw-agd_tst;
                output.Sgd_est2=S_in;
                output.adg_model2=agd_tst;
                
            end
        end

        
%% Step 6

        ff=9/lam_int;
        ff=round(ff,0);
        if ff==2
            ff=3;
        end
        
        nphy=sgolayfilt(phy,1,ff);
        
        fn={'sigma1','phi1','mu1';'sigma2','phi2','mu2';'sigma3','phi3','mu3';'sigma4','phi4','mu4';'sigma5','phi5','mu5';'sigma6','phi6','mu6';'sigma7','phi7','mu7';'sigma8','phi8','mu8';'sigma9','phi9','mu9';'sigma10','phi10','mu10';'sigma11','phi11','mu11';'sigma12','phi12','mu12';'sigma13','phi13','mu13';'sigma14','phi14','mu14';'sigma15','phi15','mu15';'sigma16','phi16','mu16'};
        
        
        drv_lam=[min(output.lam)+.5*lam_int:lam_int:max(output.lam)-.5*lam_int];
        drv2_lam=[min(output.lam)+lam_int:lam_int:max(output.lam)-lam_int];
        
        for jj=1:length(nphy)-1
            phy_drv(jj)=(nphy(jj)-nphy(jj+1))./lam_int;   %Tsai & Philpot 1998 Eq. 7
        end
        
        for jj=1:length(nphy)-2
            phy_drv2(jj)=(nphy(jj)-2*nphy(jj+1)+nphy(jj+2))/lam_int^2; %Tsai & Philpot 1998 Eq. 8
        end
        
        phy2=sgolayfilt(-phy_drv2,1,ff);
        [phi,mu,sigma]=findpeaks(phy2,drv2_lam);
        
        ind=find(sigma > 5);
        phi=phi(ind);
        mu=mu(ind);
        sigma=sigma(ind);

%% Step 7

        [Y,I]=sort(phi,'descend');
        
        [c,ia,ib]=intersect(mu,output.lam);
        
        phy_res=nphy;
        
        for ii=1:length(I)
            phi(I(ii))=phy_res(ib(I(ii)));
            gg=gauss(output.lam,sigma(I(ii)),phi(I(ii)),mu(I(ii)));
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
        
        output.number_gauss=length(phi);
        
        if length(phi) > 0
            
            for ii=1:length(phi)
                model1.model_fit.(fn{ii,1})=sigma(ii);
                model1.model_fit.(fn{ii,2})=phi(ii);
                model1.model_fit.(fn{ii,3})=mu(ii);
            end
            
            [model,model_fit,model_stats]=build_phy_model(nphy,output.lam,model1,length(phi));
            
            for ii=1:length(phi)
                model1.model_fit.(fn{ii,1})=model.(fn{ii,1});
                model1.model_fit.(fn{ii,2})=model.(fn{ii,2});
                model1.model_fit.(fn{ii,3})=model.(fn{ii,3});
            end
 
%% Step 8

            if length(phi) > 0
                if output.Sgd_est2 > 0 & output.Sgd_est2~=output.Sgd_initial_est
                    [output.final_model, output.final_model_gof, output.final_model_fitstats]=build_daisea_model(output.at_nw,output.lam,model1,length(phi),output.Sgd_est2,output.model2_lam0,output.model2_agd_lam0,start,stop);
                else
                    [output.final_model, output.final_model_gof, output.final_model_fitstats]=build_daisea_model(output.at_nw,output.lam,model1,length(phi),output.Sgd_initial_est,start,output.at_nw(ind350),start,stop);
                end
            else
                model1=[];
                if output.Sgd_est2 > 0 & output.Sgd_est2~=output.Sgd_initial_est
                    [output.final_model, output.final_model_gof, output.final_model_fitstats]=build_daisea_model(output.at_nw,output.lam,model1,length(phi),output.Sgd_est2,output.model2_lam0,output.model2_agd_lam0,start,stop);
                else
                    [output.final_model, output.final_model_gof, output.final_model_fitstats]=build_daisea_model(output.at_nw,output.lam,model1,length(phi),output.Sgd_initial_est,start,output.at_nw(ind350),start,stop);
                end
            end
            
            
            output.Sgd_final_est=output.final_model.s;
            
            yy=cdom_model_lam0_noK(output.lam,output.final_model.a,output.final_model.s,output.final_model.lam0);
            
            output.agd_final_estimate=yy;
            output.aphy_final_estimate=sgolayfilt(output.at_nw-yy,1,ff);
            
            output.final_lam0=output.final_model.lam0;
            output.agd_final_lam0=output.final_model.a;
            
                        
%% Step 00 - no output for final model            
        else
            %if we didn't satisfy the requirements, no output to signal
            %a failed attempt

            output.final_model=[];
            output.final_model_gof=[];
            output.final_model_fitstats=[];
            output.Sgd_final_est=[];
            output.agd_final_estimate=[];
            output.aphy_final_estimate=[];
            output.final_lam0=[];
            output.agd_final_lam0=[];
            output.number_gauss=[];
        end

%% Step 000 - n_points < 5, fitting not attempted
    else
        %if there are less than 5 points where the 2nd deriv = 0, do
        %not attempt to decompose the spectra
        
        output.at_nw=[];
        output.drv=[];
        output.drv2=[];
        output.n_points=[];
        output.Sgd_initial_est=[];
        output.model2_agd_lam0=[];
        output.model2_lam0=[];
        output.Sgd_est2=[];
        output.initial_adg_model=[];
        output.adg_model2=[];
        output.final_model=[];
        output.final_model_gof=[];
        output.final_model_fitstats=[];
        output.Sgd_final_est=[];
        output.agd_final_estimate=[];
        output.aphy_final_estimate=[];
        output.dg_dominant=[];
        output.final_lam0=[];
        output.agd_final_lam0=[];
        output.number_gauss=[];
    end
end

