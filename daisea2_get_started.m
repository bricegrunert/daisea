% daisea2_get_started
%
% Welcome to getting started with DAISEA2!! This script shows examples of
% how to analyze the provided example data with the daisea2 function.
% Complete documentation for daisea2 is provided at the top of that
% function. We recommend reading that documentation before proceeding
% below. 
%
% DAISEA2 requires an installation of MATLAB including the curve fitting
% and signal processing toolboxes. We also recommend that you add daisea2 
% and its subfolders to the MATLAB path. The code below assumes you have 
% done this or have the daisea2 folder as your current working directory. 
% 
% If you are unsure what toolboxes are installed, go to the MATLAB command
% window and type "ver <enter>". This will list everything that is
% installed. 
%
% We recommend running this script one section at a time, reading the
% documentaion as you go to get the best introduction to running daisea2. 
%
% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

%% Step 1 - Clean up workspace
% This will delete any variables in your workspace and close all figures.
% We recommend doing this before you run other parts of this script. 
close all
clear

%% Step 2 - Load data
% A data file titled daisea2_example_data.mat is included in the daisea2
% folder. The code below will load a data structure called "example_data"
% into your workspace. This structure includes wavelength and total
% non-water absorption (anw). 
%
% To access a field in a structure, the format is: structure.field. So, to
% access wavelength within example_data, the format would be
% example_data.wavelength
%
load('daisea2_example_data.mat');

%% Step3 - Run daisea2
% daisea2 is designed to work with one spectra at a time. You should be
% able to run it in a parallel-processing framework (parfor) if desired. We
% do not show that here. 
%
% We recommend you read the documenation at the top of the daisea2 function
% before proceeding. 
%

% This will run daisea2 with default settings (exponential CDOM/NAP model,
% negative data will be set to 0, entire spectra considered) and provide a
% figure showing each step of daisea2 fitting. 
%
% Note that we recommend, when processing large datasets, that you leave
% the doplot option to be false. Makeing a figure significantly increases
% the processing time. Also note that the figure is not closed upon the
% next iteration of each loop. Thus, if you are processing a dataset of
% 3000 spectra, this will create 3000 open figures, which usually leads to
% memory problems for MATLAB. 
for ii=1:size(example_data.anw,1)
    re(ii)=daisea2(example_data.anw(ii,:),example_data.wavelength,'doplot',true);
end

% re will be a structure array of modeled output. To access a single model
% output (corresponding to the rows of example_data.anw), for example,
% results for the 2nd spectra, the format would be re(2). 

%% This will repeat the analysis, using the hyperbolic model for CDOM/NAP
close all
for ii=1:size(example_data.anw,1)
    rh(ii)=daisea2(example_data.anw(ii,:),example_data.wavelength,'doplot',true,'cdom_model','hyperbolic');
end

%% Step 4 - plot results for both models
close all
s=get(0,'screensize');
figure;
set(gcf,'position',[round(s(3)*0.05) 1 round(s(3)*0.90) round(s(4)*0.60)],'PaperPositionMode','auto','BackingStore','off');
nr=2;
nc=5;
clrs.aph=[76 157 62]/255;
clrs.adg=[237 117 32]/255;
clrs.anw=[58 113 234]/255;

tiledlayout(nr,nc)
for ii=1:length(re)
    nexttile;
    plot(example_data.wavelength,example_data.anw(ii,:),'k-','linewidth',2);
    hold on;
    plot(re(ii).wavelength,re(ii).adg_final,'--','color',clrs.adg,'linewidth',2);
    plot(rh(ii).wavelength,rh(ii).adg_final,':','color',clrs.adg,'linewidth',2);
    plot(re(ii).wavelength,re(ii).aph_final,'--','color',clrs.aph,'linewidth',2);
    plot(rh(ii).wavelength,rh(ii).aph_final,':','color',clrs.aph,'linewidth',2);
    set(gca,'fontsize',18,'fontweight','bold');
    if rem(ii,nc)==1
        ylabel('a (m^{-1})');
    end
    if ii>nr*nc-nc
        xlabel('Wavelength (nm)');
    end
    if ii==nc
        legend({'measured a_{nw}','a_{dg} exponential','a_{dg} hyperbolic','a_{ph} exponential','a_{ph} hyperbolic'},'location','NorthEast');
        legend boxoff
    end
end