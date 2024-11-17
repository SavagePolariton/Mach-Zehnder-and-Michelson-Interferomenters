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

FONTSIZE=20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NOTICE: some lines below must be commented/uncommented depending 
%%%% on whether you are analyzing MZI or Michelson data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% dL = 27.276;  % [micron] Path length difference in the MZI2
% dL = 64.374;  % [micron] Path length difference in the MZI3
% dL = 109.609;  % [micron] Path length difference in the MZI4
% dL = 160.890;  % [micron] Path length difference in the MZI5
% dL = 224.993;  % [micron] Path length difference in the MZI6
% dL = 423.921;  % [micron] Path length difference in the MZI7
dL = 526.535;  % [micron] Path length difference in the MZI8

% dL = 12.986;  % [micron] Path length difference in the Mich2
% dL = 32.226;  % [micron] Path length difference in the Mich3
% dL = 57.879;  % [micron] Path length difference in the Mich4
% dL = 83.533;  % [micron] Path length difference in the Mich5


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load experiments
% MZI series
load("./measurements/EnricoRenzi_MZI1/09-Nov-2024 10.56.41.mat"); % calibration circuit MZI1
%
% Michelson series
% load("./measurements/EnricoRenzi_Mich1/09-Nov-2024 11.01.50.mat"); % calibration circuit Mich1
%
lambda=testResult.header.wavelength*1e-9; % data in nm
amplitude=testResult.rows.channel_2;
figure;
plot (lambda*1e6, amplitude);
title ('Calibration loopback'); 
xlabel ('Wavelength [$\mu$m]','FontSize',FONTSIZE)
ylabel ('Insertion Loss [dB]','FontSize',FONTSIZE)
hold all;

% Fit the data with a polynomial
p=polyfit((lambda-mean(lambda))*1e6, amplitude, 5);
amplitude_LOOPBACK=polyval(p,(lambda-mean(lambda))*1e6);
plot (lambda*1e6, amplitude_LOOPBACK);
% find wavelength range with usable data, in the loopback
loopback_IL = max(amplitude);
new_lambda_i=find(amplitude>loopback_IL-10);
lambda=lambda(new_lambda_i);
lambda_min = min(lambda);
lambda_max = max(lambda);
amplitude=amplitude(new_lambda_i);
% refit the loopback
LOOPBACK=polyfit((lambda-mean(lambda))*1e6, amplitude, 4);
amplitude_LOOPBACK=polyval(LOOPBACK,(lambda-mean(lambda))*1e6);
plot (lambda*1e6, [amplitude_LOOPBACK],'r-','Linewidth',5);
axis tight;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MZI data:
% load("./measurements/EnricoRenzi_MZI2/09-Nov-2024 10.55.24.mat"); % MZI4 measurement
% load("./measurements/EnricoRenzi_MZI3/09-Nov-2024 11.03.09.mat"); % MZI4 measurement
% load("./measurements/EnricoRenzi_MZI4/09-Nov-2024 10.50.16.mat"); % MZI4 measurement
% load("./measurements/EnricoRenzi_MZI5/09-Nov-2024 10.52.50.mat"); % MZI5 measurement
% load("./measurements/EnricoRenzi_MZI6/09-Nov-2024 10.59.14.mat"); % MZI6 measurement
% load("./measurements/EnricoRenzi_MZI7/09-Nov-2024 10.51.33.mat"); % MZI7 measurement
% load("./measurements/EnricoRenzi_MZI8/09-Nov-2024 11.00.31.mat"); % MZI8 measurement
%
% Michelson data
load("./measurements/EnricoRenzi_Mich2/09-Nov-2024 11.04.25.mat"); % Mich2 measurement
% load("./measurements/EnricoRenzi_Mich3/09-Nov-2024 10.54.06.mat"); % Mich3 measurement
% load("./measurements/EnricoRenzi_Mich4/09-Nov-2024 10.47.35.mat"); % Mich4 measurement
% load("./measurements/EnricoRenzi_Mich5/09-Nov-2024 10.48.57.mat"); % Mich5 measurement
lambda1=testResult.header.wavelength*1e-9; % data in nm
amplitude=testResult.rows.channel_2;
figure;
plot (lambda1*1e6, amplitude);
title ('MZI (raw data)'); 
xlabel ('Wavelength [$\mu$m]','FontSize',FONTSIZE)
ylabel ('Insertion Loss [dB]','FontSize',FONTSIZE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MZI data - calibrated
%
% data only within the bandwidth of interest.
lambda=lambda_min:min(diff(lambda1)):lambda_max;
amplitude=interp1(lambda1, amplitude, lambda,'linear');
amplitude(find(amplitude==-inf))=-50;
% calibrate data
amplitude_cal=amplitude-polyval(LOOPBACK,(lambda-mean(lambda))*1e6);

amplitude=amplitude_cal;
figure;
plot (lambda*1e6, amplitude);
title ('MZI (calibrated with loopback)'); 
xlabel ('Wavelength [$\mu$m]','FontSize',FONTSIZE)
ylabel ('Insertion Loss [dB]','FontSize',FONTSIZE)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%% uncomment this line if you work with MZI interferomenter
ng_av = mean(lambda)^2/(dL*1e-6)/fsr % for MZI

%%%%% uncomment this line if you work with Michelson interferomenter
% ng_av = mean(lambda)^2/(dL*1e-6)/fsr/2


% find starting point for curve fitting, using the ng value
% lambda0 is in microns.
lambda0 = mean(lambda) * 1e6;
n1=2.4;
%modeNumber = n1_initial * dL / lambda0 - 0.5;
%n1 = (2*floor(modeNumber)+1)*lambda0/2/dL;
n2 = (n1-ng_av)/lambda0;
n3 = -0.01;
nx_init = [n1 n2 n3];
alpha_init = 1e-3;  % propagation loss [micron^-1]
x0=[nx_init, alpha_init, 0];


% Define the MZI transfer function
% - as a Taylor expansion around the central wavelength
% - Use units of [microns] – keeps the variables closer to 1.
% - These make the curve fitting easier.
% use Matlab anonymous functions
% effective index:
neff = @(nx, lambda) ...
		(nx(1) + nx(2).*(lambda-lambda0) + nx(3).*(lambda-lambda0).^2); % MZI index

% neff([2.4, -1, 0], 1.56)  % test it.
% alpha = 1e-3;  % propagation loss [micron^-1]
% complex propagation constant
beta = @(nx, alpha, lambda) ...
		(2*pi*neff(nx, lambda)./lambda - 1i*alpha/2*ones(1,length(lambda)) );
% beta([2.4, -1, 0], 1e-3, [1.56, 1.57]) % test it.

%%%%% uncomment this line if you work with MZI interferomenter
T_MZI = @(X, lambda) ...
        (10*log10( 0.25* abs(1+exp(-1i*beta(X(1:3), X(4), lambda)*dL)).^2) +X(5) ); % MZI Transfer Function

%%%%% uncomment this line if you work with Michelson interferomenter
% T_MZI = @(X, lambda) ...
        % (10*log10( 0.5* abs(1+exp(-1i*beta(X(1:3), X(4), lambda)*dL*2)).^2) +X(5) ); % Michelson Transfer Function



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

%%%%% uncomment this line if you work with NZI interferomenter
n_shift = lambda_offset*lambda0/fsr/dL; % NZI
%%%%% uncomment this line if you work with Michelson interferomenter
% n_shift = lambda_offset*lambda0/fsr/dL/2; % Mich

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


% legend
% legend('(470,223.1)nm','(470,215.3)nm','(615,223.1)nm','(615,215.3)nm','Mich3 LB','Mich4 LB','Mich5 LB','Mich3 BL','Mich4 BL','Mich5 BL','Nominal')
% legend('(470,223.1)nm','(470,215.3)nm','MZI2','MZI3','MZI4','MZI5','MZI6','MZI7','MZI8','(520,223.1)nm','(520,215.3)nm','Nominal')