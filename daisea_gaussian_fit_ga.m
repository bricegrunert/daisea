function [res,cnt,yhat]=daisea_gaussian_fit_ga(absorption,wavelength,ngauss,gaussian_fit,lambda_start,lambda_stop)

% [model, model_gof, model_fitstats]=build_daisea_model(absorption,wavelength,lambda_start,lambda_stop,lam0,model,ngauss,peak_lim,Sdg_est,alam0)
%
% Function identifies a suitable model to fit total non-water absorption
% spectra using a single exponential and pre-determined number of Gaussian 
% curves. Steps to initially estimate input values derived from steps in
% daisea.m
%
% Inputs:
%     absorption   = input absorption spectra (or other signal) used to
%                    optimize input Gaussian components from a least
%                    squares approach, format = vector
%     wavelength   = wavelength values affiliated with absorption spectra,
%                    format = vector
%     model        = model object with initial Gaussian estimates
%     ngauss       = number of Gaussian components used in model
%     Sdg_est      = input estimate of NAP/CDOM absorption spectral slope
%     lam0         = wavelength used to initialize fitting of NAP/CDOM
%                    absorption
%     alam0        = input estimate of NAP/CDOM absorption at lam0
%     lambda_start = minimum wavelength considered (e.g. 350 nm),
%                    must be within spectral range of wavelength input
%     lambda_stop  = maximum wavelength considered (e.g. 700 nm),
%                    must be within spectral range of wavelength input
%
% Returns:
%     model          = output from Gaussian decomposition of
%                      absorption
%     model_gof      = goodness of fit for model
%     model_fitstats = fitting stats for model
% 

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

%% Housekeeping
if nargin==4
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

% find peak locations
pl=NaN(1,ngauss);
for ii=1:ngauss
    [jumk,ind]=min(abs(x-gaussian_fit.(['mu' num2str(ii)])));
    pl(ii)=ind;
end
clear ind

%% generate first set of parents
if true % parents generated randomly between limits
    parents=gaussian_fit;
    fn={'sigma','phi','mu'};
    for ii=1:ngauss
        for ff=1:length(fn)
            parents.([fn{ff} num2str(ii)])=[parents.([fn{ff} num2str(ii)]); lb.([fn{ff} num2str(ii)]) + (ub.([fn{ff} num2str(ii)])-lb.([fn{ff} num2str(ii)])).*rand(g.population-1,1)];
        end
    end
elseif false % parents generated as +/- of initial guess
    lpercent=50;
    upercent=150;
    parents=gaussian_fit;
    fn={'sigma','phi','mu'};
    for ii=1:ngauss
        for ff=1:length(fn)
            if lb.([fn{ff} num2str(ii)])==ub.([fn{ff} num2str(ii)])
                parents.([fn{ff} num2str(ii)])=repmat(parents.([fn{ff} num2str(ii)]),g.population,1);
            else
                asdf=(lpercent + (upercent-lpercent).*rand(g.population-1,1))./100;
                parents.([fn{ff} num2str(ii)])=[parents.([fn{ff} num2str(ii)]); repmat(parents.([fn{ff} num2str(ii)]),g.population-1,1).*asdf];
            end
        end
    end
else
    % generate half random and half close to first guess
    parents=gaussian_fit;
    fn={'sigma','phi','mu'};
    for ii=1:ngauss
        for ff=1:length(fn)
            parents.([fn{ff} num2str(ii)])=[parents.([fn{ff} num2str(ii)]); lb.([fn{ff} num2str(ii)]) + (ub.([fn{ff} num2str(ii)])-lb.([fn{ff} num2str(ii)])).*rand(g.population./2-1,1)];
        end
    end
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
end

% run parents through function
lt690=wavelength<=690;
yhat=daisea_gaussian_model_ga(ngauss,x,parents);
% estimate error
tst=sum((repmat(y,g.population,1)-yhat).^2,2);
%tst2=sum((repmat(y,g.population,1)-yhat)>0,2)./size(yhat,2);
tst3=sum((repmat(y(pl),g.population,1)-yhat(:,pl)).^2,2);
%tst4=sum((repmat(y(pl(end)),g.population,1)-yhat(:,pl(end))).^2,2);
%tst5=sum((repmat(y,g.population,1)-yhat_cdom)<=0,2)==0;
%tst5=sum((repmat(y(:,lt690),g.population,1)-yhat_cdom(:,lt690))>0,2)./size(yhat_cdom(:,lt690),2)==1;
tst6=sum(yhat(:,lt690)>0,2)./size(yhat(:,lt690),2)==1;

if false
    close all;
    figure;
    set(gcf,'position',[584   406   798   595],'PaperPositionMode','auto','BackingStore','off');%,'PaperOrientation','landscape');
    plot(x,yhat,'-');
    hold on;
    plot(x,y,'k-','linewidth',4);
    plot(x,yhat(1,:),'b-','linewidth',4);
    %plot(x,yhat_cdom,'y-');
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
    %r=((1./tst)./max(1./tst)) + (tst2./max(tst2)) + ((1./tst3)./max(1./tst3));
    %r=((1./tst)./max(1./tst)) + ((1./tst3)./max(1./tst3)) + ((1./tst4)./max(1./tst4));
    r=(((1./tst)./max(1./tst)) + ((1./tst3)./max(1./tst3))) + tst6;
    r=r./max(r);
    %r=(1./tst3)./max(1./tst3);
    %r=(tst+tst3);
    %r=(1/r)./(max(1/r));
    tournament=[randsample(g.population,g.population/2,true,r) randsample(g.population,g.population/2,true,r) randsample(g.population,g.population/2,true,r) randsample(g.population,g.population/2,true,r)];
    %tournament=[randsample(g.population,g.population/2) randsample(g.population,g.population/2) randsample(g.population,g.population/2) randsample(g.population,g.population/2)];
    [junk,ind]=sort(r(tournament),2);
    useme=NaN(g.population/2,2);
    for ii=1:g.population/2
        useme(ii,:)=tournament(ii,ind(ii,end-1:end));
    end
    for ff=1:length(fn)
        children.(fn{ff})=mean(parents.(fn{ff})(useme),2);
    end
    if false
        yhat=daisea_gaussian_model_ga(ngauss,x,children);
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
        yhat=daisea_gaussian_model_ga(ngauss,x,children);
        plot(x,yhat(end-g.population/2:end,:),'r-','linewidth',1);
    end    
    % reassess fitness and cull population
    yhat=daisea_gaussian_model_ga(ngauss,x,children);
    tst=sum((repmat(y,g.population*2,1)-yhat).^2,2);
    %tst2=sum((repmat(y,g.population*2,1)-yhat)>0,2)./size(yhat,2);
    tst3=sum((repmat(y(pl),g.population*2,1)-yhat(:,pl)).^2,2);
    %tst4=sum((repmat(y(pl(end)),g.population*2,1)-yhat(:,pl(end))).^2,2);
    %tst5=sum((repmat(y,g.population*2,1)-yhat_cdom)<=0,2)==0;
    %tst5=sum((repmat(y(:,lt690),g.population*2,1)-yhat_cdom(:,lt690))>0,2)./size(yhat_cdom(:,lt690),2)==1;
    tst6=sum(yhat(:,lt690)>0,2)./size(yhat(:,lt690),2)==1;
    r=((1./tst)./max(1./tst)) + ((1./tst3)./max(1./tst3)) + tst6;
    r=r./max(r);
    [junk,ind]=sort(r,'descend');
    %[junk,ind]=sort(tst);
    for ff=1:length(fn)
        parents.(fn{ff})=children.(fn{ff})(ind(1:g.population));
    end
    yhat=yhat(ind(1:g.population),:);
    tst=tst(ind(1:g.population));
    %tst2=tst2(ind(1:g.population));
    tst3=tst3(ind(1:g.population));
    %tst4=tst4(ind(1:g.population));
    %tst5=tst5(ind(1:g.population));
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
