function [res,cnt,yhat]=daisea_adg_ga(absorption,wavelength,cdom_model,cdom_fit,percent_dg,wlind,lambda_start,lambda_stop,offset)

% [output,qcflag]=daisea_adg_ga(output,wlind)
%
% Genetic algorithm to refine initial estimate of adg. 
%
% Inputs:
%  output     = output structrue from daisea
%  wlind      = index structure for wavelength positions
%  cdom_model = exponential, exponential_k, stretched_exponential,
%               hyperbolic
%  cdom_fit   = fit object from daisea_initial_adg (includes a and s)
%  percent_dg = percent contribution of adg at 440 nm
%  wlind      = structure of indexed locations for specific wavelengths
%  lambda_start = lowest wavelength in fit
%  lambda_stop  = highest wavelength in fit
%  offset       = data offset (for samples <10% aph at 440 nm)
%

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

% Created on 2022-09-22 by AC

%% Housekeeping
if nargin==6
    lambda_start=min(wavelength);
    lambda_stop=max(wavelength);
end

ind=wavelength>=lambda_start & wavelength<=lambda_stop;

absorption=absorption(:);
wavelength=wavelength(:);

x=wavelength(ind)';
y=absorption(ind)';

%% genetic algorithm choices and parameter settings
clear g
g.population=50;  % number of coefficient sets in population, must be even
g.generations=50; % number of generations to grow
g.mutation=100;     % mutation rate, percent, must be a whole number

% upper and lower bounds for parameters
clear lb ub
lb.percent_dg=percent_dg-0.1;
ub.percent_dg=percent_dg+0.1;
if lb.percent_dg<0
    lb.percent_dg=0.01;
end
if ub.percent_dg>1
    ub.percent_dg=0.99;
end

[junk,useme]=min(abs(x-440));
switch cdom_model
    case 'exponential'
        %lind=find(x==cdom_fit.lam0);
        lb.a=0;
        ub.a=y(useme);
        lb.s=0;
        ub.s=0.03;
        lb.lam0=cdom_fit.lam0;
        ub.lam0=cdom_fit.lam0;
    case 'hyperbolic'
        lb.a=0;
        %lb.a=cdom_fit.a*0.1;
        %ub.a=cdom_fit.a*1.3;
        ub.a=y(useme);
        clear junk useme
        lb.s=0;
        ub.s=cdom_fit.s+4;
        %lb.a=cdom_fit.a*0.2;
        %ub.a=cdom_fit.a*1.2;
        %lb.s=cdom_fit.s-1;
        %ub.s=cdom_fit.s+1;
end
ub.offset=offset;
lb.offset=offset;

%% generate first set of parents
fn=fieldnames(cdom_fit);
for ff=1:length(fn)
    if strcmp(fn{ff},'lam0')
        parents.(fn{ff})=repmat(cdom_fit.(fn{ff}),g.population,1);
    else
        parents.(fn{ff})=[cdom_fit.(fn{ff}); lb.(fn{ff}) + (ub.(fn{ff})-lb.(fn{ff})).*rand(g.population-1,1)];
        %parents.(fn{ff})=[repmat(cdom_fit.(fn{ff}),g.population./2+1,1); lb.(fn{ff}) + (ub.(fn{ff})-lb.(fn{ff})).*rand(g.population./2-1,1)];
    end
end
parents.percent_dg=[percent_dg; lb.percent_dg + (ub.percent_dg-lb.percent_dg).*rand(g.population-1,1)];
parents.offset=repmat(offset,size(parents.percent_dg));

% run parents through function
yhat=daisea_combined_model_ga(cdom_model,0,x,parents);
yhatphy=repmat(y,g.population,1)-yhat;

% estimate error
lt600=wavelength<=600;
lt690=wavelength<=650;
if any(y(lt690)<0) || any(y(lt690)==0)
    lt690=wavelength<=500;
end
% Minimize sum of squares difference, be as close to anw as possible
tst=sum((repmat(y(:,lt600),g.population,1)-yhat(:,lt600)).^2,2);
% Want estimated curve all below anw (ie. aph is positive)
tst2=sum(yhatphy(:,lt690)>0,2)./size(yhat(:,lt690),2)==1;
% Want Ratio of aph(350)/aph(440) is <1.5
if isnan(wlind.ind350)
    tst3=yhatphy(:,wlind.ind400)./yhatphy(:,wlind.ind440)<=1.2;
else
    tst3=yhatphy(:,wlind.ind350)./yhatphy(:,wlind.ind440)<=1.5;
end
% Want estimated adg above 0 (ie. no negative absorption)
tst4=sum(yhat(:,lt690)>0,2)/size(yhat(:,lt690),2)==1;

if false
    close all;
    figure;
    set(gcf,'position',[584   406   798   595],'PaperPositionMode','auto','BackingStore','off');%,'PaperOrientation','landscape');
    plot(x,yhat,'-');
    hold on;
    plot(x,y,'k-','linewidth',4);
    plot(x,yhat(1,:),'b-','linewidth',4);
    set(gca,'fontsize',18,'fontweight','bold');
    xlabel('Wavelength (nm)')
    ylabel('a_{nw} (m^{-1})')
end

%% genetic algorithm
cnt=0;
fn=fieldnames(parents);
while cnt<g.generations
    cnt=cnt+1;
    % assess fitness and generate children
    clear children
    r=(((1./tst)./max(1./tst))./1 + tst3 + tst4) .* tst2;
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
        yhat=daisea_combined_model_ga(cdom_model,0,x,children);
        plot(x,yhat,'g-','linewidth',1);
        drawnow
    end
    % mutate
    chance=randi(100,g.population,1);
    chance=chance<=g.mutation;
    for ff=1:length(fn)
        p.(fn{ff})=(80 + (120-80).*rand(g.population,1))./100;
        p.(fn{ff})(~chance)=1;
    end
    useme=randsample(g.population,g.population/2,false);
    for ff=1:length(fn)
        children.(fn{ff})=[parents.(fn{ff}); children.(fn{ff}); parents.(fn{ff})(useme).*p.(fn{ff})(useme)];
        children.(fn{ff})(children.(fn{ff})<lb.(fn{ff}))=lb.(fn{ff});
        children.(fn{ff})(children.(fn{ff})>ub.(fn{ff}))=ub.(fn{ff});
    end
    if false
        yhat=daisea_combined_model_ga(cdom_model,0,x,children);
        plot(x,yhat(end-g.population/2:end,:),'r-','linewidth',1);
        drawnow
    end    
    % reassess fitness and cull population
    yhat=daisea_combined_model_ga(cdom_model,0,x,children);
    yhatphy=repmat(y,g.population*2,1)-yhat;
    tst=sum((repmat(y(:,lt600),g.population*2,1)-yhat(:,lt600)).^2,2);
    tst2=sum(yhatphy(:,lt690)>0,2)./size(yhat(:,lt690),2)==1;
    if isnan(wlind.ind350)
        tst3=yhatphy(:,wlind.ind400)./yhatphy(:,wlind.ind440)<=1.2;
    else
        tst3=yhatphy(:,wlind.ind350)./yhatphy(:,wlind.ind440)<=1.5;
    end
    tst4=sum(yhat(:,lt690)>0,2)/size(yhat(:,lt690),2)==1;
    r=(((1./tst)./max(1./tst)) + tst3 + tst4) .* tst2;
    r=r./max(r);
    [junk,ind]=sort(r,'descend');
    for ff=1:length(fn)
        parents.(fn{ff})=children.(fn{ff})(ind(1:g.population));
    end
    yhat=yhat(ind(1:g.population),:);
    tst=tst(ind(1:g.population));
    tst2=tst2(ind(1:g.population));
    tst3=tst3(ind(1:g.population));
    tst4=tst4(ind(1:g.population));
    clear r
end

%%
if false
    %close all;
    figure;
    set(gcf,'position',[584   406   798   595],'PaperPositionMode','auto','BackingStore','off');%,'PaperOrientation','landscape');
    plot(x,yhat,'-');
    hold on;
    plot(x,y,'k-','linewidth',4);
    plot(x,yhat(1,:),'r-','linewidth',4);
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
