function y=daisea_combined_model_ga(cdom_model,ngauss,x,c)

% y=daisea_combined_model_ga(cdom_model,ngauss,x,c)
%
% Apply DAISEA models for combined fits with genetic algorithm. Returns a
% vector of predicted values. If ngauss=0, only CDOM/NAP estimate will be
% returned. 
%
% Inputs:
%  cdom_model = 'exponential' or 'hyperbolic'
%  ngauss     = number of gaussian components, if set to 0 will only return
%               CDOM results
%  x          = wavelength
%  c          = structure of fitted parameters (a, s, lam0, mu1, phi1,
%               sigma1, etc...). If ngauss=0, then only a, s and lam0 are
%               used. 
% 

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

% Created on 2023-03-20 by AC

%%
% If not specified, percent_dg is set to 1
if ~isfield(c,'percent_dg')
    c.percent_dg=1;
end
% If not specified, offeset is set to 0
if ~isfield(c,'offset')
    c.offset=0;
end

% Retrieve CDOM
switch cdom_model
    case 'exponential'
        y=(c.percent_dg.*c.a.*exp(-(x-c.lam0).*c.s))-c.offset;
    case 'hyperbolic'
        %y=(c.percent_dg.*c.a.*(x./532).^-c.s)-c.offset;
        y=(c.percent_dg.*c.a.*(x./440).^-c.s)-c.offset;
    case 'stretched_exponential'
        y=(c.percent_dg.*c.a.*exp(-(c.s.*(x-c.lam0)).^c.beta));
    case 'exponential_k'
        y=(c.percent_dg.*c.a.*exp(-(x-c.lam0).*c.s)+c.k);
end

% Retrieve Gaussian
switch ngauss
    case 1
        y=y+c.phi1.*exp(-(x-c.c.mu1).^2./(2.*c.c.sigma1.^2));
    case 2
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2));
    case 3
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2));
    case 4
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2));
    case 5
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2));
    case 6
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2));
    case 7
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2));
    case 8
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2));
    case 9
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2));
    case 10
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2));
    case 11
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2));
    case 12
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2));
    case 13
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2))+c.phi13.*exp(-(x-c.mu13).^2./(2.*c.sigma13.^2));
    case 14
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2))+c.phi13.*exp(-(x-c.mu13).^2./(2.*c.sigma13.^2))+c.phi14.*exp(-(x-c.mu14).^2./(2.*c.sigma14.^2));
    case 15
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2))+c.phi13.*exp(-(x-c.mu13).^2./(2.*c.sigma13.^2))+c.phi14.*exp(-(x-c.mu14).^2./(2.*c.sigma14.^2))+c.phi15.*exp(-(x-c.mu15).^2./(2.*c.sigma15.^2));
    case 16
        y=y+c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2))+c.phi13.*exp(-(x-c.mu13).^2./(2.*c.sigma13.^2))+c.phi14.*exp(-(x-c.mu14).^2./(2.*c.sigma14.^2))+c.phi15.*exp(-(x-c.mu15).^2./(2.*c.sigma15.^2))+c.phi16.*exp(-(x-c.mu16).^2./(2.*c.sigma16.^2));
end