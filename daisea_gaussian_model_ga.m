function y=daisea_gaussian_model_ga(ngauss,x,c)

% DAISEA models for gaussian fits with genetic algorithm.
%

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

% Created on 2024-11-25 by AC

%%
% Retrieve Gaussian
switch ngauss
    case 1
        y=c.phi1.*exp(-(x-c.c.mu1).^2./(2.*c.c.sigma1.^2));
    case 2
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2));
    case 3
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2));
    case 4
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2));
    case 5
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2));
    case 6
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2));
    case 7
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2));
    case 8
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2));
    case 9
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2));
    case 10
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2));
    case 11
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2));
    case 12
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2));
    case 13
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2))+c.phi13.*exp(-(x-c.mu13).^2./(2.*c.sigma13.^2));
    case 14
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2))+c.phi13.*exp(-(x-c.mu13).^2./(2.*c.sigma13.^2))+c.phi14.*exp(-(x-c.mu14).^2./(2.*c.sigma14.^2));
    case 15
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2))+c.phi13.*exp(-(x-c.mu13).^2./(2.*c.sigma13.^2))+c.phi14.*exp(-(x-c.mu14).^2./(2.*c.sigma14.^2))+c.phi15.*exp(-(x-c.mu15).^2./(2.*c.sigma15.^2));
    case 16
        y=c.phi1.*exp(-(x-c.mu1).^2./(2.*c.sigma1.^2))+c.phi2.*exp(-(x-c.mu2).^2./(2.*c.sigma2.^2))+c.phi3.*exp(-(x-c.mu3).^2./(2.*c.sigma3.^2))+c.phi4.*exp(-(x-c.mu4).^2./(2.*c.sigma4.^2))+c.phi5.*exp(-(x-c.mu5).^2./(2.*c.sigma5.^2))+c.phi6.*exp(-(x-c.mu6).^2./(2.*c.sigma6.^2))+c.phi7.*exp(-(x-c.mu7).^2./(2.*c.sigma7.^2))+c.phi8.*exp(-(x-c.mu8).^2./(2.*c.sigma8.^2))+c.phi9.*exp(-(x-c.mu9).^2./(2.*c.sigma9.^2))+c.phi10.*exp(-(x-c.mu10).^2./(2.*c.sigma10.^2))+c.phi11.*exp(-(x-c.mu11).^2./(2.*c.sigma11.^2))+c.phi12.*exp(-(x-c.mu12).^2./(2.*c.sigma12.^2))+c.phi13.*exp(-(x-c.mu13).^2./(2.*c.sigma13.^2))+c.phi14.*exp(-(x-c.mu14).^2./(2.*c.sigma14.^2))+c.phi15.*exp(-(x-c.mu15).^2./(2.*c.sigma15.^2))+c.phi16.*exp(-(x-c.mu16).^2./(2.*c.sigma16.^2));
end