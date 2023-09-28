% This script is to fit the initital, exponential growth phase of model
% pathogen load to LCMV viral load. This estimates the initial serum
% pathogen load that fits the data, P(t=0), for both acute and chronic
% pathogens, as explained in the SI Text
%
% Date of creation: October 16th, 2022
%
% (C) Sam Jamaleddine 2022

%% Get everything ready

% Load data
load('LCMVserumWherry.mat')

time_acute = serum.acute.mod(:,1);
time_chronic = serum.chronic.mod(:,1);
virus_acute = serum.acute.mod(:,2);
virus_chronic = serum.chronic.mod(:,2);

% Use the data from the exponential growth phase to estimate betas/ICs
ta = time_acute(1:3); % there's not enough data in acute viral load
Va = virus_acute(1:3);

tc = time_chronic(1:5);
Vc = virus_chronic(1:5);

% Fit
linFit_acute = fit(ta,log(Va),'poly1');
linFit_chronic = fit(tc,log(Vc),'poly1');

% Plot
pars_acute = coeffvalues(linFit_acute);
pars_chronic = coeffvalues(linFit_chronic);

figurer(1)
subplotter(1,2,1)
plot(linFit_acute,ta,log(Va))

subplotter(1,2,2)
plot(linFit_chronic,tc,log(Vc))

beta_a = pars_acute(1); % estimate for acute replication rate
V0_a = exp(pars_acute(2)); % estimate for initial LCMV-Arm pathogen load
beta_c = pars_chronic(1); % estimate for chronic replication rate
V0_c = exp(pars_chronic(2)); % estimate for initial LCMV-Cl13 load
