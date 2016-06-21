function [power,freq,maxp]=periodog(x,fontsize)
%PERIODOG Periodogram of a time series.
%   [POWER,FREQ,MAXP]=PERIODOG(X) returns the POWER and FREQuency such that 
%   PLOT(FREQ,POWER) is the periodogram (estimate of the spectral density). 
%   Period (i.e. 1/FREQ) MAXP of the maximum value of the periodogram is 
%   also returned.
%   If no output parameters are supplied the periodogram is plotted 
%   automatically. PERIODOG(X,FONTSIZE) allows to specify a fontsize 
%   different than 14.
%
%   References:
%   [1] P.Stoica, R.L.Moses (1997) Introduction to Spectral Analysis,
%   Prentice Hall.
%   [2] R.Weron (2006) Modeling and Forecasting Electricity Loads and 
%   Prices: A Statistical Approach, Wiley.


%   Written by Rafal Weron (2000.12.27)
%   Revised by Rafal Weron (2006.09.23)

if nargin<2, 
    fontsize = 14; 
end;

if length(x)/2 ~= floor(length(x)/2), 
    x = x(1:end-1);
end;

% Compute the Fourier transform of the data
Y = fft(x);
N = length(Y);
% Remove the first component of Y, which is simply the sum of the data
Y(1) = [];
% Compute the power as the squared absolute value of Y
% A plot of power versus frequency is the 'periodogram'
power = (abs(Y(1:N/2)).^2)/length(x);
% Define the frequencies
nyquist = 1/2;
freq = ((1:N/2)/(N/2)*nyquist)';
period = 1./freq;
% Find the maximum value of the periodogram
[max_p ind_p] = max(power);
maxp = period(ind_p);

if nargout<1,
    % Plot the periodogram
    subplot(2,1,1)
    plot(freq,power,'.-')
    grid on
    set(gca,'Box','on','fontsize',fontsize);
    xlabel('Frequency','fontsize',fontsize);
    ylabel('Power','fontsize',fontsize);
    % A plot of power vs. 1/freq indicates the length of cycles
    subplot(2,1,2)
    plot(1./freq,power,'.-')
    grid on
    set(gca,'Box','on','fontsize',fontsize);
    xlabel('Period','fontsize',fontsize);
    ylabel('Power','fontsize',fontsize);
end;
