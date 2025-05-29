function output=daisea2(absorption,wavelength,varargin)

% daisea2    Decompose total non-water absorption into CDOM/NAP and
%            phytoplankton components. 
%
% output = daisea2(absorption,wavelength) returns an output structure
% containing fit results and confidence intervals. 
%
% DAISEA2 takes total non-water absorption and decomposes it into CDOM/NAP
% and phytoplankton components using a series of genetic algorithms.
% Gaussian decomposition of estimated phytoplankton absorption allows for
% the assessment of pigment composition. 
%
% DAISEA2 is intended for hyperspectral data from 350-700nm with band 
% width intervals of 1nm to 5nm. Updates from the original DAISEA now allow
% for analysis of spectra starting at or near 400 nm (e.g. ACS data). Data 
% at or near wavelegnths of 350 (or 400, if 350 is missing), 440, 555, 680 
% and 690 nm are specifically used. The nearest match to these wavelengths 
% is utilized if an exact match is not found. 
%
% DAISEA2 requires installation of the curve fitting and signal processing
% toolboxes. 
%
% Associated manuscript:
% Grunert, B.K., A.B. Ciochetto and C.B. Mouw (2025) A hyperspectral 
% approach for retrieving inherent optical properties, 
% phytoplankton pigments, and associated uncertainties from non-water 
% absorption. ******** INSERT DOI AND REST OF CITATION *********
%
% Please refer to this paper for a detailed description of each step with 
% an accompanying flow chart as well as providing a thorough description of
% model output.
%
% Data requirements:
%   - DAISEA2 requires that all data have a consistent wavelength
%     interval for the entire spectra. This should be a whole number
%     between 1nm and 5nm. If they do not (e.g. ACS data), then, by default
%     the function will interpolate the data to a wavelength interval of
%     floor(median(diff(wavelength))). Other interpolation options are 
%     available via the 'wavelength_interval', name/value input argument 
%     described below. 
%
%   - DAISEA2 requires all data be positive. Options for how negative data
%     are treated can be set via the 'negative_data' name/value input
%     argrument described below. 
%
%
% Inputs:
%  absorption = total non-water absorption in m-1 (vector)
%  wavelength = wavelength in nm (vector)
%
%
% Optional name/value pairs:
%  cdom_model          = 'exponential' or 'hyperbolic'. Default=exponential
%  negative data       = DAISEA2 requires all data to be positive. Negative
%                        data can either be 'offset' (add lowest value to
%                        entire dataset to make it positive), 'zero' (set
%                        to 0), or 'remove' (set to NaN and ignored),
%                        default is 'zero'.
%  wavelength_interval = DAISEA2 requires all data to have a regular
%                        wavelength interval. If they do not (e.g. ACS
%                        data), then, by default the program will
%                        interpolate the data to a wavelength interval of
%                        floor(median(diff(wavelength))). Otherwise, the
%                        user may enter a value here to override the
%                        default.
%
%                        If your data have a regular wavelength
%                        interval and you enter a value here that is
%                        different from that interval, then the data will
%                        be interpolated to the interval specified here.
%
%                        If your data have a regular wavelength
%                        interval but get removed (e.g. because they are
%                        negative or NaN), then the truncated data set
%                        will be interpolated back to the original
%                        wavelength interval. Again, if you enter a value
%                        here that is different from the wavelength
%                        interval of your data, then the data will be
%                        interpolated to the interval specified here.
%  minimum_wavelength  = minimum wavelength included in fit. Data below 
%                        this value will be ignored. 
%                        Default=min(wavelength)
%  maximum_wavelength  = maximum wavelength inclued in fit. Data above 
%                        this value will be ignored. 
%                        Default=max(wavelength)
%  doplot              = true/false, default=false. Show figure of fitting
%                        steps
%  notify              = true/false, default=false. Print warnings to 
%                        screen.
%
%
% Returns:
%     output = structure containing information from derivative analysis
%     and iterative spectral fitting or from full DAISEA2 algorithm. If a
%     field is unspecified, it is left blank in the structure. Fields are
%     as follows:
%           wavelength           = input wavelength (nm)
%           anw                  = input total non-water absorption (m-1)
%           percent_ph           = estimated contibution of phytoplankton 
%                                  at 440 nm
%           percent_dg           = estimated contribution of CDOM and NAP 
%                                  at 440 nm
%           wavelength_interval  = wavelength spread between data points (nm)
%           min_wavelength       = minimum wavelength included in fit (nm)
%           max_wavelength       = maximum wavelength included in fit (nm)
%           cdom_model           = CDOM model used for fit of adg
%           drv                  = first derivative of anw
%           drv2                 = second derivative of anw
%           n_points             = number of data points where drv2 = 0
%           adg_fit1             = first estimate of adg, fitted parameters
%           adg_model1           = vector of first estimate of adg
%           adg_model1_adjusted  = adg_model1 adjusted for percent_dg
%           offset               = for samples with percent_ph<10%, an
%                                  offest is applied to shift data above
%                                  690nm to 0. For all samples with
%                                  percent_ph>10%, offset is 0. 
%           adg_fit2             = second estimate of adg via genetic 
%                                  algorithm, fitted parameters
%           adg_model2           = vector of second estimate of adg
%           aph_step5            = vector of estimated aph after adg_model2
%                                  (aph_step5 = anw - adg_model2)
%           number_gauss         = number of Gaussian components used to 
%                                  parameterize aph in final_model
%           generations          = number of generations in final fit
%           ga_reps              = fitted values for each GA itteration
%           ga_fits              = vector of anw for each entry in ga_reps
%           final_model          = final fitted parameters
%           final_model_ci       = confidence intervals for fitted
%                                  parameters
%           adg_final            = final estimate of adg
%           aph_final            = anw - adg_final
%
%
% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

%% Step 0
% Re-assign input absorption and wavelength to appear in output structure in
% the proper format (single column for each parameter), remove NaN's and
% negative values, assign start and stop lambda if needed

% Create settings structure and check on user input
p=inputParser;
addParameter(p,'notify',false,@(x) islogical(x));
addParameter(p,'minimum_wavelength',min(wavelength),@(x) isnumeric(x) && isscalar(x) && (x>0));
addParameter(p,'maximum_wavelength',max(wavelength),@(x) isnumeric(x) && isscalar(x) && (x>0));
addParameter(p,'negative_data','zero',@(x) any(validatestring(x,{'offset','zero','remove'})));
addParameter(p,'wavelength_interval',[],@(x) isnumeric(x) && isscalar(x) && (x>0));
addParameter(p,'cdom_model','exponential',@(x) any(validatestring(x,{'exponential','exponential_k','stretched_exponential','hyperbolic'})));
addParameter(p,'doplot',false,@(x) islogical(x));
parse(p,varargin{:});
settings=p.Results;

% Make sure absorption and wavelength inputs are not flipped
if any(absorption>=350) && sum(wavelength<350)==length(wavelength)
    error('ERROR DAISEA2: Wavelength is all under 350 nm. Check order of absorption and wavelength inputs.')
end

% Ensure proper vector orientation
% Initialize return structure
output=struct('wavelength',wavelength(:),'anw',absorption(:),'percent_ph',[],'percent_dg',[],...
    'wavelength_interval',[],'min_wavelength',[],'max_wavelength',[],...
    'cdom_model',settings.cdom_model,'drv',[],'drv2',[],'n_points',[],...
    'adg_fit1',[],'adg_model1',[],'adg_model1_adjusted',[],'offset',0,....
    'adg_fit2',[],'adg_model2',[],'aph_step5',[],...
    'number_gauss',[],'generations',[],'ga_reps',[],'ga_fits',[],...
    'final_model',[],'final_model_ci',[],...
    'adg_final',[],...
    'aph_final',[],...
    'mu1',[],'phi1',[],'sigma1',[],'mu2',[],'phi2',[],'sigma2',[]);

% initialize figure
if settings.doplot
    s=get(0,'screensize');
    figure;
    set(gcf,'position',[round(s(3)*0.1) 1 round(s(3)*0.80) round(s(4)*0.66)],'PaperPositionMode','auto','BackingStore','off');%,'PaperOrientation','landscape');
    fcnt=0;
    nr=2;
    nc=3;
end

% Remove NaNs, take care of negative values, set wavelength interval
output=daisea_wavelength_interval_check(output,settings);

% Find data falling within the desired fitting window (e.g. 350-700 nm)
[~,start_ind]=min(abs(output.wavelength-settings.minimum_wavelength));
[~,stop_ind]=min(abs(output.wavelength-settings.maximum_wavelength));

output.min_wavelength=output.wavelength(start_ind);
output.max_wavelength=output.wavelength(stop_ind);

output.wavelength=output.wavelength(start_ind:stop_ind);
output.anw=output.anw(start_ind:stop_ind);

% Look for specific wavelength locaitons
wlind=daisea_find_wavelengths(output.wavelength,settings.notify);

% empirical relationship from training dataset - Grunert et al.
% submitted June 2018

% Retrieve band ratio 555/680
ratio=output.anw(wlind.ind555)/output.anw(wlind.ind680);
if ratio > 0.685
    output.percent_ph=1.038*exp(-.9257*(output.anw(wlind.ind555)/output.anw(wlind.ind680)));
else
    output.percent_ph=2.088*exp(-1.946*(output.anw(wlind.ind555)/output.anw(wlind.ind680)));
end
if output.percent_ph>1
    output.percent_ph=0.5;
end
output.percent_dg=1-output.percent_ph;

%% move smoothing to begining
if true
    ff=9/output.wavelength_interval;
    ff=round(ff,0);
    if ff==2
        ff=3;
    end
    output.anw=sgolayfilt(output.anw,1,ff);
end

%% Steps 1 and 2
% Get initial estimate of adg
[output,qcflag,ind]=daisea_initial_adg(output,wlind,settings.cdom_model);
if qcflag; return; end

if settings.doplot
    fcnt=fcnt+1;
    subplot(nr,nc,fcnt)
    plot(output.wavelength,output.anw,'k-','linewidth',1);
    hold on;
    plot(output.wavelength(ind),output.anw(ind),'b.','markersize',10);
    plot(output.wavelength,output.adg_model1,'r--','linewidth',1);
    plot(output.wavelength,output.adg_model1_adjusted,'r:','linewidth',1);
    set(gca,'fontsize',18,'fontweight','bold','xlim',[min(output.wavelength) max(output.wavelength)]);
    xlabel('Wavelength (nm)');
    ylabel('a_{nw} or a_{dg} (m^{-1})');
    clear asdf
    ylim=get(gca,'ylim');
    title('Steps 1 & 2');
    legend({'a_{nw}','Points used in fit','Initial model','Model adjusted for %a_{ph}(440)'},'location','NorthEast');
    legend boxoff
end

%% Samples below 10%
if true
    if output.percent_ph<0.1
        % are there points where estimated adg is above measured at_nw?
        wltstmax=690;
        tst=output.adg_model1_adjusted(output.wavelength<=wltstmax)-output.anw(output.wavelength<=wltstmax);
        if any(tst>0)
            output.offset=max(tst);
            %output.anw=output.anw+max(tst);
        end
    end
end

%% Steps 3 - 5
nreps=10;
clear res n yhat res2 res3
n=NaN(1,nreps);
yhat=NaN(nreps,length(ind));
for ii=1:nreps
    [res(ii),n(ii),yhat(ii,:)]=daisea_adg_ga(output.anw(ind),output.wavelength(ind),output.cdom_model,output.adg_fit1,output.percent_dg,wlind,output.min_wavelength,output.max_wavelength,output.offset);
end
fn=fieldnames(res);
tst=cell2mat(squeeze(struct2cell(res)));
for ii=1:length(fn)
    res2.(fn{ii})=tst(ii,:)';
    res3.(fn{ii})=median(res2.(fn{ii}));
end
output.adg_fit2=rmfield(res3,'percent_dg');
output.adg_model2=daisea_combined_model_ga(output.cdom_model,0,output.wavelength,res3);
output.aph_step5=output.anw-output.adg_model2;

output.percent_dg=res3.percent_dg;
output.percent_ph=1-output.percent_dg;

if settings.doplot
    fcnt=fcnt+1;
    subplot(nr,nc,fcnt)
    plot(output.wavelength,output.anw,'k-','linewidth',1);
    hold on;
    plot(output.wavelength,output.adg_model2,'r-','linewidth',1);
    plot(output.wavelength,output.aph_step5,'g-','linewidth',1);
    set(gca,'fontsize',18,'fontweight','bold','xlim',[min(output.wavelength) max(output.wavelength)],'ylim',ylim);
    xlabel('Wavelength (nm)');
    ylabel('a (m^{-1})');
    title('Steps 3 - 5');
    legend({'a_{nw}','a_{dg}','a_{ph}'},'location','NorthEast');
    legend boxoff
    drawnow
end

%% Step 6
res=daisea_gauss_peak_locations(output.wavelength,output.aph_step5);
ff=res.ff;
nphy=res.nphy;
number_gauss=res.number_gauss;
sigma=res.sigma;
mu=res.mu;
phi=res.phi;

output.mu1=res.mu1;
output.mu2=mu;
output.phi1=res.phi1;
output.phi2=phi;
output.sigma1=res.sigma1;
output.sigma2=sigma;

clear res

fn={'sigma','phi','mu'};

% if we didn't satisfy the requirements, no output to signal
% a failed attempt
if isempty(phi)
    return;
end

output.number_gauss=number_gauss;
r.sigma=sigma;
r.mu=mu;
r.phi=phi;

if settings.doplot
    fcnt=fcnt+1;
    subplot(nr,nc,fcnt)
    plot(output.wavelength,output.aph_step5,'g-','linewidth',1);
    hold on;
    for ii=1:length(r.mu)
        [junk,p]=min(abs(output.wavelength-r.mu(ii)));
        plot(r.mu(ii),output.aph_step5(p),'b.','markersize',12);
        clear junk p
    end
    set(gca,'fontsize',18,'fontweight','bold','xlim',[min(output.wavelength) max(output.wavelength)]);%,'ylim',ylim);
    xlabel('Wavelength (nm)');
    ylabel('a_{ph} (m^{-1})');
    title('Steps 6 and 7');
    legend({'a_{ph}','Peak locations'},'location','NorthEast');
    legend boxoff
    drawnow
end

% do fit
clear gaussian_fit
for ii=1:output.number_gauss
    for jj=1:3
        gaussian_fit.([fn{jj} num2str(ii)])=r.(fn{jj})(ii);
    end
end
clear r

[model,model_gof,model_fitstats]=daisea_gaussian_fit(output.aph_step5,output.wavelength,gaussian_fit,output.number_gauss);

for ii=1:output.number_gauss
    for jj=1:3
        gaussian_fit.([fn{jj} num2str(ii)])=model.([fn{jj} num2str(ii)]);
    end
end

%% Step 8
nreps=10;
clear res n yhat res2 res3
n=NaN(1,nreps);
yhat=NaN(nreps,length(output.wavelength));
for ii=1:nreps
    [res(ii),n(ii),yhat(ii,:)]=daisea_combined_fit_ga(output.anw,output.wavelength,output.number_gauss,gaussian_fit,output.cdom_model,output.adg_fit2,output.percent_dg,output.min_wavelength,output.max_wavelength);
end
fn=fieldnames(res);
tst=cell2mat(squeeze(struct2cell(res)));
for ii=1:length(fn)
    res2.(fn{ii})=tst(ii,:)';
    res3.(fn{ii})=median(res2.(fn{ii}));
    res4.(fn{ii})=[max(res2.(fn{ii})) min(res2.(fn{ii}))];
end
output.generations=n;
output.ga_reps=res2;
output.ga_fits=yhat;
output.final_model=res3;
output.final_model_ci=res4;

output.adg_final=daisea_combined_model_ga(output.cdom_model,0,output.wavelength,output.final_model);
output.aph_final=sgolayfilt(output.anw-output.adg_final,1,ff);

if settings.doplot
    %figure(1)
    fcnt=fcnt+1;
    subplot(nr,nc,fcnt)
    plot(output.wavelength,output.anw,'k-','linewidth',1);
    hold on;
    plot(output.wavelength,output.adg_final,'r-','linewidth',1);
    plot(output.wavelength,output.aph_final,'g-','linewidth',1);
    set(gca,'fontsize',18,'fontweight','bold','xlim',[min(output.wavelength) max(output.wavelength)],'ylim',ylim);
    xlabel('Wavelength (nm)');
    ylabel('a (m^{-1})');
    title('Step 8');
    legend({'a_{nw}','a_{dg}','a_{ph}'},'location','NorthEast');
    legend boxoff
    drawnow
    pause(1)
end

