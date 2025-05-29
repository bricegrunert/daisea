function yhat = daisea_apply_cdom_model(wm,x,varargin)

% yhat = daisea_apply_cdom_model(wm,x,varargin)
%
% Retreive yhat based on initial fitting results. 
%
% Inputs:
%  wm = which model, 'exponential' or 'hyperbolic'
%  x  = wavelength vector
%
% Name/argument pairs:
%  fit = fit object
%  a = intercept
%  s = slope
%  beta = parameter for stretched exponential
%  lam0 = lamda 0
%  k = offset for exponential_k
%  percent_dg = percent contribution of adg at 440 nm
%

% Copyright (c) 2025 Brice K. Grunert, Audrey B. Ciochetto, Colleen B. Mouw
% See license.txt
% email: b.grunert@csuohio.edu, a.ciochetto@csuohio.edu, cmouw@uri.edu

%%
p=inputParser;
addParameter(p,'fit',[],@(x) isobject(x) || isstruct(x));
addParameter(p,'a',[],@(x) isnumeric(x) && isscalar(x));
addParameter(p,'s',[],@(x) isnumeric(x) && isscalar(x));
addParameter(p,'beta',[],@(x) isnumeric(x) && isscalar(x));
addParameter(p,'lam0',[],@(x) isnumeric(x) && isscalar(x));
addParameter(p,'k',[],@(x) isnumeric(x) && isscalar(x));
addParameter(p,'percent_dg',1,@(x) isnumeric(x) && isscalar(x));
parse(p,varargin{:});
p=p.Results;

switch wm 
    case 'exponential'
        if ~isempty(p.fit)
            yhat=cdom_model_exponential(x,p.fit.a.*p.percent_dg,p.fit.s,p.fit.lam0);
        else
            yhat=cdom_model_exponential(x,p.a.*p.percent_dg,p.s,p.lam0);
        end
    case 'exponential_k'
        if ~isempty(p.fit)
            yhat=cdom_model_exponential_k(x,p.fit.a.*p.percent_dg,p.fit.s,p.fit.lam0,p.fit.k);
        else
            yhat=cdom_model_exponential_k(x,p.a.*p.percent_dg,p.s,p.lam0,p.k);
        end
    case 'stretched_exponential'
        if ~isempty(p.fit)
            yhat=cdom_model_stretched_exponential(x,p.fit.a.*p.percent_dg,p.fit.s,p.fit.beta,p.fit.lam0);
        else
            yhat=cdom_model_stretched_exponential(x,p.a.*p.percent_dg,p.s,p.beta,p.lam0);
        end
    case 'hyperbolic'
        if ~isempty(p.fit)
            yhat=cdom_model_hyperbolic(x,p.fit.a.*p.percent_dg,p.fit.s);
        else
            yhat=cdom_model_hyperbolic(x,p.a.*p.percent_dg,p.s);
        end
end