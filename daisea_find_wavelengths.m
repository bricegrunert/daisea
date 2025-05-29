function ind=daisea_find_wavelengths(wavelength,notif)

% ind=daisea_find_wavelengths(wavelength,notif)
% 
% Find closest wavelengths to those used throughout DAISEA processing. 
%

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

% Created on 2022-09-19 by AC

%%

[junk,ind.ind555]=min(abs(wavelength-555));
if junk>5
    ind.ind555=NaN;
    if notif
        fprintf('\nWARNING DAISEA: Closest wavelength to 555 nm is %0.2f nm away.',junk);
    end
end

[junk,ind.ind680]=min(abs(wavelength-680));
if junk>5
    ind.ind680=NaN;
    if notif
        fprintf('\nWARNING DAISEA: Closest wavelength to 680 nm is %0.2f nm away.',junk);
    end
end

[junk,ind.ind690]=min(abs(wavelength-690));
if junk>5
    ind.ind690=NaN;
    if notif
        fprintf('\nWARNING DAISEA: Closest wavelength to 690 nm is %0.2f nm away.',junk);
    end
end

[junk,ind.ind350]=min(abs(wavelength-350));
if junk>5
    ind.ind350=NaN;
    if notif
        fprintf('\nWARNING DAISEA: Closest wavelength to 350 nm is %0.2f nm away.',junk);
    end
end

[junk,ind.ind440]=min(abs(wavelength-440));
if junk>5
    ind.ind440=NaN;
    if notif
        fprintf('\nWARNING DAISEA: Closest wavelength to 440 nm is %0.2f nm away.',junk);
    end
end

[junk,ind.ind400]=min(abs(wavelength-400));
if junk>5
    ind.ind400=NaN;
    if notif
        fprintf('\nWARNING DAISEA: Closest wavelength to 400 nm is %0.2f nm away.',junk);
    end
end