function [h,conf_lo,conf_hi]=gph_weron(x,cutoff,conf_level,fontsize)
%GPH Geweke-Porter-Hudak estimator of the Hurst exponent.
%   H=GPH(X) returns the Hurst exponent H of a time series X estimated
%   using the Geweke-Porter-Hudak (1983) spectral estimator for periods 

%   lower than max(period)^CUTOFF, where CUTOFF=0.5.
%   H=GPH(X,CUTOFF) allows to specify a CUTOFF different then 0.5.
%   [H,CONF_LO,CONF_HI]=GPH(X,CUTOFF,CONF_LEVEL) also returns  
%   the asymptotic confidence interval [CONF_LO,CONF_HI] at a given 
%   (two sided) confidence level (default: CONF_LEVEL=0.95).
%
%   If there are no output parameters, log(power) is automatically 
%   plotted against log(4*sin(period/2)^2). GPH(...,CONF_LEVEL,FONTSIZE) 
%   allows to specify a fontsize different than 14 in the plotted figure.
%
%   References:
%   [1] J.Geweke, S.Porter-Hudak (1983) The estimation and application of 
%   long memory time series models, Journal of Time Series Analysis 4, 
%   221-238.
%   [2] R.Weron (2002) Estimating long range dependence: finite sample 
%   properties and confidence intervals, Physica A 312, 285-299.

%	Written by Rafal Weron (2000.12.27)
%	Revised by Rafal Weron (2001.02.04,2011.09.30)

if nargin<4, 
    fontsize = 14; 
end;
if nargin<3, 
    conf_level = 0.95;
end;
if nargin<2, 
    cutoff = 0.5;
end;

x = x(:);
[power,freq] = periodog(x);
N = length(x);
T = find(freq<N^-cutoff);

% Compute the Hurst exponent
PT = polyfit(log(4*sin(freq(T)/2).^2),log(power(T)),1);
h = 0.5 - PT(1);

% Compute confidence intervals
sig = pi/(sqrt(6*(length(T)-1)) * std(log(4*sin(freq(T)/2).^2)));
conf_lo = norminv((1-conf_level)/2,.5,sig);
conf_hi = norminv(1-(1-conf_level)/2,.5,sig);

% Plot power vs. frequency on a loglog scale
if nargout<1,
    plot(log10(4*sin(freq(T)/2).^2),log10(power(T)),'.');
    hold on
    maxx = max((4*sin(freq(T)/2).^2));
    minx = min((4*sin(freq(T)/2).^2));
    plot(log10([maxx minx]),log10(([maxx minx].^PT(1))*exp(PT(2))),'r')
    hold off

    set(gca,'Box','on','fontsize',fontsize);
    xlabel('log_{10} (4 sin^2(freq/2))','fontsize',fontsize);
    ylabel('log_{10} (power)','fontsize',fontsize);
    disp(['Hurst exponent: ' num2str(h)])
end;