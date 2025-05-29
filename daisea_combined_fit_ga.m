function [res,cnt,yhat]=daisea_combined_fit_ga(absorption,wavelength,ngauss,gaussian_fit,cdom_model,cdom_fit,percent_dg,lambda_start,lambda_stop)

% [res,cnt,yhat]=daisea_combined_fit_ga(absorption,wavelength,ngauss,gaussian_fit,cdom_model,cdom_fit,percent_dg,lambda_start,lambda_stop)
%
% Genetic algorithm for a combined fit of CDOM/NAP and gaussian components.
% 
% Inputs:
%     absorption   = anw (vector)
%     wavelength   = wavelength (vector)
%     ngauss       = number of gaussian components
%     gaussian_fit = fit object with initial fitted values for each peak
%                    (mu, sigma, phi)
%     cdom_model   = 'exponential' or 'hyperbolic'
%     cdom_fit     = fit object of initial parameters for CDOM/NAP fit
%     percent_dg   = percent contribution of adg at 440 nm
%     lambda_start = minimum wavelength considered (e.g. 350 nm),
%                    must be within spectral range of wavelength input
%     lambda_stop  = maximum wavelength considered (e.g. 700 nm),
%                    must be within spectral range of wavelength input
%
% Returns:
%     res  = structure of fitted parameters
%     cnt  = number of generations
%     yhat = predicted anw
% 

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

%% Housekeeping
if nargin==7
    lambda_start=min(wavelength);
    lambda_stop=max(wavelength);
end

ind=wavelength>=lambda_start & wavelength<=lambda_stop;

absorption=absorption(:);
wavelength=wavelength(:);

x=wavelength(ind)';
y=absorption(ind)';

%% genetic algorithm choices and parameter settings
% fit mu or fix it at locaitons provided?
fit_mu=false;
fit_pdg=false;

clear g
g.population=200;  % number of coefficient sets in population, must be even
g.generations=200; % number of generations to grow
g.mutation=100;     % mutation rate, percent, must be a whole number

% Try changing how percent_dg is applied in final fit
if false % Results are slighly better the old way. 2023-09-11 AC
    cdom_fit.a=cdom_fit.a.*percent_dg;
    percent_dg=1;
end

% upper and lower bounds for parameters
clear lb ub
fn=fieldnames(gaussian_fit);
% get good upper bounds for phi
for ii=1:ngauss
    ub.(['phi' num2str(ii)])=interp1(x,y,gaussian_fit.(['mu' num2str(ii)]));
end
for ff=1:length(fn)
    switch fn{ff}(1)
        case 's' % sigma
            lb.(fn{ff})=5;
            ub.(fn{ff})=50;
        case 'p' %phi
            lb.(fn{ff})=0;
            %ub.(fn{ff})=gaussian_fit.(fn{ff})*2;
        case 'm' % mu
            if fit_mu
                lb.(fn{ff})=gaussian_fit.(fn{ff})-10;
                ub.(fn{ff})=gaussian_fit.(fn{ff})+10;
            else
                lb.(fn{ff})=gaussian_fit.(fn{ff});
                ub.(fn{ff})=gaussian_fit.(fn{ff});
            end
    end
end
switch cdom_model
    case 'exponential'
        [junk,lind]=min(abs(x-cdom_fit.lam0));
        lb.a=0;
        ub.a=y(lind);
        lb.s=cdom_fit.s-0.002;
        ub.s=cdom_fit.s+0.003;
        lb.lam0=cdom_fit.lam0;
        ub.lam0=cdom_fit.lam0;
        lb.offset=cdom_fit.offset;
        ub.offset=cdom_fit.offset;
    case 'hyperbolic'
        [junk,lind]=min(abs(x-440));
        lb.a=0;
        %ub.a=cdom_fit.a*1.3;
        ub.a=y(lind);
        lb.s=cdom_fit.s-2;
        if lb.s<0; lb.s=0; end
        ub.s=cdom_fit.s+5;
        lb.offset=cdom_fit.offset;
        ub.offset=cdom_fit.offset;
    case 'exponential_k'
        lind=find(x==cdom_fit.lam0);
        lb.a=0;
        ub.a=y(lind);
        lb.s=cdom_fit.s-0.002;
        ub.s=cdom_fit.s+0.003;
        g.k_lb=0;
        g.l_ub=1;
        lb.lam0=cdom_fit.lam0;
        ub.lam0=cdom_fit.lam0;
    case 'stretched_exponential'
        lind=find(x==cdom_fit.lam0);
        lb.a=0;
        ub.a=y(lind);
        lb.s=cdom_fit.s-0.002;
        ub.s=cdom_fit.s+0.003;
        lb.lam0=cdom_fit.lam0;
        ub.lam0=cdom_fit.lam0;
end
if fit_pdg
    lb.percent_dg=percent_dg-0.1;
    ub.percent_dg=percent_dg+0.1;
    if lb.percent_dg<0; lb.percent_dg=0.01; end
    if ub.percent_dg>1; lb.percent_dg=0.99; end
else
    lb.percent_dg=percent_dg;
    ub.percent_dg=percent_dg;
end

% find peak locations
pl=NaN(1,ngauss);
for ii=1:ngauss
    [jumk,ind]=min(abs(x-gaussian_fit.(['mu' num2str(ii)])));
    pl(ii)=ind;
end
clear ind

%% generate first set of parents

% generate half random and half close to first guess
parents=gaussian_fit;
fn={'sigma','phi','mu'};
for ii=1:ngauss
    for ff=1:length(fn)
        parents.([fn{ff} num2str(ii)])=[parents.([fn{ff} num2str(ii)]); lb.([fn{ff} num2str(ii)]) + (ub.([fn{ff} num2str(ii)])-lb.([fn{ff} num2str(ii)])).*rand(g.population./2-1,1)];
    end
end
fn=fieldnames(cdom_fit);
for ff=1:length(fn)
    if strcmp(fn{ff},'lam0')
        parents.(fn{ff})=repmat(cdom_fit.(fn{ff}),g.population./2,1);
    else
        parents.(fn{ff})=[cdom_fit.(fn{ff}); lb.(fn{ff}) + (ub.(fn{ff})-lb.(fn{ff})).*rand(g.population./2-1,1)];
    end
end
parents.percent_dg=[percent_dg; lb.percent_dg + (ub.percent_dg-lb.percent_dg).*rand(g.population./2-1,1)];
% half close
lpercent=50;
upercent=150;
fn={'sigma','phi','mu'};
for ii=1:ngauss
    for ff=1:length(fn)
        if lb.([fn{ff} num2str(ii)])==ub.([fn{ff} num2str(ii)])
            parents.([fn{ff} num2str(ii)])=[parents.([fn{ff} num2str(ii)]); repmat(gaussian_fit.([fn{ff} num2str(ii)]),g.population./2,1)];
        else
            asdf=(lpercent + (upercent-lpercent).*rand(g.population./2,1))./100;
            parents.([fn{ff} num2str(ii)])=[parents.([fn{ff} num2str(ii)]); repmat(gaussian_fit.([fn{ff} num2str(ii)]),g.population./2,1).*asdf];
        end
    end
end
fn=fieldnames(cdom_fit);
for ff=1:length(fn)
    if lb.(fn{ff})==ub.(fn{ff})
        parents.(fn{ff})=[parents.(fn{ff}); repmat(cdom_fit.(fn{ff}),g.population./2,1)];
    else
        %parents.(fn{ff})=[cdom_fit.(fn{ff}); lb.(fn{ff}) + (ub.(fn{ff})-lb.(fn{ff})).*rand(g.population-1,1)];
        asdf=(lpercent + (upercent-lpercent).*rand(g.population./2,1))./100;
        parents.(fn{ff})=[parents.(fn{ff}); repmat(cdom_fit.(fn{ff}),g.population./2,1).*asdf];
    end
end
parents.percent_dg=[parents.percent_dg; lb.percent_dg + (ub.percent_dg-lb.percent_dg).*rand(g.population./2,1)];

% run parents through function
lt690=wavelength<=690;
yhat=daisea_combined_model_ga(cdom_model,ngauss,x,parents);
yhat_cdom=daisea_combined_model_ga(cdom_model,0,x,parents);
% estimate error
tst=sum((repmat(y,g.population,1)-yhat).^2,2);
tst3=sum((repmat(y(pl),g.population,1)-yhat(:,pl)).^2,2);
tst5=sum((repmat(y(:,lt690),g.population,1)-yhat_cdom(:,lt690))>0,2)./size(yhat_cdom(:,lt690),2)==1;
tst6=sum(yhat_cdom(:,lt690)>0,2)./size(yhat_cdom(:,lt690),2)==1;

if false
    close all;
    figure;
    set(gcf,'position',[584   406   798   595],'PaperPositionMode','auto','BackingStore','off');%,'PaperOrientation','landscape');
    plot(x,yhat,'-');
    hold on;
    plot(x,y,'k-','linewidth',4);
    plot(x,yhat(1,:),'b-','linewidth',4);
    plot(x,yhat_cdom,'y-');
    set(gca,'fontsize',18,'fontweight','bold');
    xlabel('Wavelength (nm)')
    ylabel('a_{nw} (m^{-1})')
    drawnow
end

%% genetic algorithm
cnt=0;
fn=fieldnames(parents);
while (tst(1)>0.0005 || cnt<100) && cnt<g.generations
    cnt=cnt+1;
    % assess fitness and generate children
    clear children
    r=(((1./tst)./max(1./tst)) + ((1./tst3)./max(1./tst3))) + tst5 + tst6;
    r=r./max(r);
    tournament=[randsample(g.population,g.population/2,true,r) randsample(g.population,g.population/2,true,r) randsample(g.population,g.population/2,true,r) randsample(g.population,g.population/2,true,r)];
    [junk,ind]=sort(r(tournament),2);
    useme=NaN(g.population/2,2);
    for ii=1:g.population/2
        useme(ii,:)=tournament(ii,ind(ii,end-1:end));
    end
    for ff=1:length(fn)
        children.(fn{ff})=mean(parents.(fn{ff})(useme),2);
    end
    if false
        yhat=daisea_combined_model_ga(cdom_model,ngauss,x,children);
        plot(x,yhat,'g-','linewidth',1);
    end
    % mutate
    chance=randi(100,g.population,1);
    chance=chance<=g.mutation;
    for ff=1:length(fn)
        p.(fn{ff})=(90 + (110-90).*rand(g.population,1))./100;
        p.(fn{ff})(~chance)=1;
    end
    useme=randsample(g.population,g.population/2,false);
    for ff=1:length(fn)
        children.(fn{ff})=[parents.(fn{ff}); children.(fn{ff}); parents.(fn{ff})(useme).*p.(fn{ff})(useme)];
        children.(fn{ff})(children.(fn{ff})<lb.(fn{ff}))=lb.(fn{ff});
        children.(fn{ff})(children.(fn{ff})>ub.(fn{ff}))=ub.(fn{ff});
    end
    if false
        yhat=daisea_combined_model_ga(cdom_model,ngauss,x,children);
        plot(x,yhat(end-g.population/2:end,:),'r-','linewidth',1);
    end    
    % reassess fitness and cull population
    yhat=daisea_combined_model_ga(cdom_model,ngauss,x,children);
    yhat_cdom=daisea_combined_model_ga(cdom_model,0,x,children);
    tst=sum((repmat(y,g.population*2,1)-yhat).^2,2);
    tst3=sum((repmat(y(pl),g.population*2,1)-yhat(:,pl)).^2,2);
    tst5=sum((repmat(y(:,lt690),g.population*2,1)-yhat_cdom(:,lt690))>0,2)./size(yhat_cdom(:,lt690),2)==1;
    tst6=sum(yhat_cdom(:,lt690)>0,2)./size(yhat_cdom(:,lt690),2)==1;
    r=((1./tst)./max(1./tst)) + ((1./tst3)./max(1./tst3)) + tst5 + tst6;
    r=r./max(r);
    [junk,ind]=sort(r,'descend');
    for ff=1:length(fn)
        parents.(fn{ff})=children.(fn{ff})(ind(1:g.population));
    end
    yhat=yhat(ind(1:g.population),:);
    yhat_cdom=yhat_cdom(ind(1:g.population),:);
    tst=tst(ind(1:g.population));
    tst3=tst3(ind(1:g.population));
    tst5=tst5(ind(1:g.population));
    tst6=tst6(ind(1:g.population));
    clear r
end

%%
if false
    %close all;
    figure;
    set(gcf,'position',[584   406   798   595],'PaperPositionMode','auto','BackingStore','off');%,'PaperOrientation','landscape');
    plot(x,yhat,'-');
    hold on;
    plot(x,y,'k-','linewidth',2);
    plot(x,yhat(1,:),'r-','linewidth',2);
    set(gca,'fontsize',18,'fontweight','bold');
    xlabel('Wavelength (nm)')
    ylabel('a_{nw} (m^{-1})')
    drawnow
end

% get best result
for ff=1:length(fn)
    res.(fn{ff})=parents.(fn{ff})(1);
end
if false
    close all;
    figure;
    plot(x,yhat,'-');
    hold on;
    plot(x,yhat(ind(1),:),'r-',x,y,'k-','linewidth',2);
else
    %plot(x,yhat(ind(1),:),'r-','linewidth',4);
end

% see parents
if false
    lastgen=[];
    for ff=1:length(fn)
        lastgen(:,ff)=parents.(fn{ff});
    end
end

%%
yhat=yhat(1,:);