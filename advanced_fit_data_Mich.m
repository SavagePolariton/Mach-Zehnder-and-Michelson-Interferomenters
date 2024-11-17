% Phot1x_fit_wg_compactmodel.m
% by Lukas Chrostowski, 2015
% adapted by Enrico Maria Renzi, 2024

clc
clear all
% close all
%
%
set(groot,'defaulttextinterpreter','latex');
set(0,'defaultAxesFontSize',16)
set(0,'defaultLegendFontSize',16)
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

POLY_ORDER=4;

% dL = 12.986;  % [micron] Path length difference in the Mich2
% dL = 32.226;  % [micron] Path length difference in the Mich3
% dL = 57.879;  % [micron] Path length difference in the Mich4
dL = 83.533;  % [micron] Path length difference in the Mich5


% Michelson series
% load("./measurements/EnricoRenzi_Mich2/09-Nov-2024 11.04.25.mat"); % Mich2 measurement
% load("./measurements/EnricoRenzi_Mich3/09-Nov-2024 10.54.06.mat"); % Mich3 measurement
% load("./measurements/EnricoRenzi_Mich4/09-Nov-2024 10.47.35.mat"); % Mich4 measurement
load("./measurements/EnricoRenzi_Mich5/09-Nov-2024 10.48.57.mat"); % Mich5 measurement

lambda=testResult.header.wavelength*1e-9; % data in nm
amplitude=testResult.rows.channel_2;

% Curve fit data to a polynomial for baseline correction
p=polyfit((lambda-mean(lambda))*1e6, amplitude, 4);
amplitude_baseline=polyval(p,(lambda-mean(lambda))*1e6); 

% Perform baseline correction to flatten the spectrum
% Use the curve polynomial, and subtract from original data
amplitude_corrected = amplitude - amplitude_baseline;
amplitude_corrected = amplitude_corrected + max(amplitude_baseline) - max(amplitude);

% data only within the wavelength range of interest.
lambda_min = min(lambda);	% Can limit the analysis to a range of wavelengths
lambda_max = max(lambda);   %  if the data on the edges is noisy
lambda_min = 1.54e-6;
% lambda_max = 1.57e-6;
lambda1=lambda_min:min(diff(lambda)):lambda_max;
amplitude=interp1(lambda, amplitude_corrected, lambda1,'linear');
lambda=lambda1;
amplitude(find(amplitude==-inf))=-30;  % check if there are -infinity data points

% plot baseline corrected spectrum
figure;
plot (lambda*1e6, amplitude);
xlabel ('Wavelength [$\mu$m]');
ylabel ('Transmission [dB]');
axis tight
title ('Experimental data (baseline corrected, wavelength range)');

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

% Find ng from autocorrelation-based frequency estimation of spectrum
% auto-correction
[r,lags]=xcorr(amplitude); 
r=r(ge(lags,0));
lags=lags(ge(lags,0));
figure
plot(lags,r);
% estimate the frequency
d=diff(r);
start = find(gt(d,0)); start=start(1);
[peak_m, peak_i]=max(r(start:end));
peak_i=peak_i+start;  % location of the 1st peak in the autocorrelation
hold on;
plot(peak_i,0,'s','MarkerSize',20);
title ('Autocorrelation of spectrum')
xlabel('lag, sample number');

fsr = peak_i * mean(diff(lambda))
ng_av = mean(lambda)^2/(dL*1e-6)/fsr/2 % note the 2 factor!


% find starting point for curve fitting, using the ng value
% lambda0 is in microns.
lambda0 = mean(lambda) * 1e6;
n1=2.4;
%modeNumber = n1_initial * dL / lambda0 - 0.5;
%n1 = (2*floor(modeNumber)+1)*lambda0/2/dL;
n2 = (n1-ng_av)/lambda0;
nx_init = [n1 n2 0];
alpha_init = 1e-3;  % propagation loss [micron^-1]
x0=[nx_init, alpha_init, 0];


% Define the MZI transfer function
% - as a Taylor expansion around the central wavelength
% - Use units of [microns] – keeps the variables closer to 1.
% - These make the curve fitting easier.
% use Matlab anonymous functions
% effective index:
neff = @(nx, lambda) ...
		(nx(1) + nx(2).*(lambda-lambda0) + nx(3).*(lambda-lambda0).^2); 
% neff([2.4, -1, 0], 1.56)  % test it.
% alpha = 1e-3;  % propagation loss [micron^-1]
% complex propagation constant
beta = @(nx, alpha, lambda) ...
		(2*pi*neff(nx, lambda)./lambda - 1i*alpha/2*ones(1,length(lambda)) );
% beta([2.4, -1, 0], 1e-3, [1.56, 1.57]) % test it.
% MZI transfer function
T_MZI = @(X, lambda) ...
        (10*log10( 0.5* abs(1+exp(-1i*beta(X(1:3), X(4), lambda)*2*dL)).^2) +X(5) );
% T_MZI([2.4, -1, 0, 1e-3], [1.56, 1.57]) % test it.


figure;
plot (lambda*1e6, amplitude);
hold all;
plot (lambda0, -40,'s','MarkerSize',20);
plot(lambda*1e6, T_MZI(x0, lambda*1e6),'--','LineWidth',3);
xlabel ('Wavelength [$\mu$m]');
ylabel ('Transmission [dB]');
axis tight
title ('MZI model (initial parameters)');


% Autocorrelation again, to find the shift between the fit function and experimental data
[r,lags]=xcorr(amplitude, T_MZI(x0, lambda*1e6)); 
r=r(ge(lags,0));
lags=lags(ge(lags,0));
[peak_m, peak_i]=max(r);
lambda_offset = peak_i(1) * mean(diff(lambda));
n_shift = lambda_offset*lambda0/fsr/dL/2;
x0(1)=x0(1)+n_shift;

figure;
plot (lambda*1e6, amplitude);
hold all;
plot (lambda0, -40,'s','MarkerSize',20);
plot(lambda*1e6, T_MZI(x0, lambda*1e6),'--','LineWidth',3);
xlabel ('Wavelength [$\mu$m]');
ylabel ('Transmission [dB]');
axis tight
title ('MZI model (initial parameters, with shift)');


% Curve fit:  
[xfit,resnorm] = lsqcurvefit(T_MZI,x0,lambda*1e6,amplitude);
xfit
r=corrcoef(amplitude,T_MZI(xfit, lambda*1e6));
r2=r(1,2).^2

figure;
plot (lambda*1e6, amplitude);
hold all;
plot(lambda*1e6, T_MZI(xfit, lambda*1e6),'LineWidth',3); 
xlabel ('Wavelength [$\mu$m]');
ylabel ('Transmission [dB]');
axis tight
title ('MZI model (fit parameters)');

% Check if the fit is good.  If so, find ng
if (ge(r2,0.8))
  % plot ng curve
  figure;
  hold on
  neff_fit = neff(xfit(1:3),lambda*1e6);
  dndlambda=diff(neff_fit)./diff(lambda); dndlambda=[dndlambda, dndlambda(end)];
  ng=(neff_fit - lambda .* dndlambda);
  plot(lambda*1e6, ng, 'LineWidth',4);
  xlabel ('Wavelength [$\mu$m]');
  ylabel ('Group index, $n_g$');
  axis tight
  title ('Group index (from MZI fit)');
    
  % waveguide parameters at lambda0
  ng0 = xfit(1) - lambda0*xfit(2)
end