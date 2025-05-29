function output=daisea_wavelength_interval_check(output,settings)

% output = daisea_wavelength_interval_check(output,settings)
%
% Data must be NaN and non-negative for DAISEA. After removing NaNs and
% modifying negative values, this function does a wavelength interval
% check. If data are removed during this step, thereby creating an
% inconsistent wavelength interval, remaining data will be interpolated
% back to the original wavelength interval. So, if your data were 1 nm
% apart, but some were NaN or negative and were removed, then these holes
% will be filled in so the final data product is again 1 nm apart with no
% gaps. If data were removed at the min/max wavelengths of the dataset,
% these are not interpolated. Instead the dataset will be trimmed to the
% new min/max. 
%
% DAISEA requires a wavelength intervals <= 5 nm. User can specify 
% wavelength interval in settings structure. This function checks the 
% wavelength interval and interpolates data as per user input or with the 5
% nm restriction. 
%
% Note that the wavelength_interval field of the settings structure
% superceeds all other calculations/requirements. So, if a value is
% specified in the structure, data will be interpolated to that wavelength
% interval regardless of a lack of data gaps. 
%
% Inputs:
%  output     = output structure from daisea. Uses wavelength and anw.
%  settings   = settings structure as defined in daisea function. Used here
%               are:
%               negative_data = specify how to modify negative values
%               (offset, zero, remove), see daisea function. 
%               wavelength_interval = allows user to specify wavelength
%               interval. Data are interpolated using interp1 with a linear
%               interpolation.
%               notify = print to screen what's going on (true/false)
% 

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

% Created on 2022-09-19 by AC
 
%%
% Find wavelength interval from raw data
output.wavelength_interval=unique(diff(output.wavelength));

% Get rid of NaNs
output.wavelength=output.wavelength(~isnan(output.anw));
output.anw=output.anw(~isnan(output.anw));
if isempty(output.anw)
    error('ERROR DAISEA: Data are all NaN.');
end

% If minimum is less than 0, add offset or set to zero
if min(output.anw)<0
    switch settings.negative_data
        case 'offset'
            output.anw=output.anw-min(output.anw);
        case 'zero'
            output.anw(output.anw<0)=0;
        case 'remove'
            output.wavelength=output.wavelength(output.anw>0);
            output.anw=output.anw(output.anw>0);
    end
end

% If wavelength interval is not consistent, fix
lam_int2=unique(diff(output.wavelength));
fixme=false;
if (length(output.wavelength_interval)==1 && length(lam_int2)>1)
    fixme=true; % fix but use original output.wavelength interval
elseif length(output.wavelength_interval)>1 && length(lam_int2)>1
    fixme=true; % fix and determine output.wavelength interval
    output.wavelength_interval=floor(median(diff(output.wavelength)));
end
% get wavelength interval
if any((~isnan(settings.wavelength_interval) | ~isempty(settings.wavelength_interval)) & settings.wavelength_interval~=output.wavelength_interval)
    fixme=true;
    output.wavelength_interval=settings.wavelength_interval;
end
if output.wavelength_interval>5
    fprintf('\nWARNING DAISEA: wavelength interval is %0.2f nm. Data will be interpolated to the maximum allowable interval of 5 nm.',output.wavelength_interval);
    fixme=true;
    output.wavelength_interval=5;
end
% Put onto new grid
if fixme
    if settings.notify
        fprintf('\nDAISEA: Data are being interpolated to a %0.2f nm wavelength interval',output.wavelength_interval);
    end
    wavelength2=ceil(min(output.wavelength)):output.wavelength_interval:floor(max(output.wavelength));
    absorption2=interp1(output.wavelength,output.anw,wavelength2,'linear',NaN);
    output.wavelength=wavelength2(:);
    output.anw=absorption2(:);
else
    output.wavelength=output.wavelength;
    output.anw=output.anw;
end
clear lam_int2 fixme wavelength2 absorption2