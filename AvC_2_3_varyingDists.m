%% This script determines effect of reactivity mode or span (Figs. 3, S6)
% This script takes some time to run- the results were therefore saved as
% .mat files, and loaded in AvC_figs_main or AvC_figs_supp to plot the
% results
%
%   AvC_2_3_effectOfDists.mat = effect of mode and span plotted in Fig 3 of
%   the main manuscript
%
%   AvC_2_3_effectOfMode_uniform = effect of varying mode, when exhaustion
%   rates are instead drawn from uniform distributions (Fig. S6D). The span
%   was not saved but can check that this is also the same.
% 
% (C) Hassan (Sam) Jamaleddine, 2022

%% Do the computations

n = 1; % number of simulations at each bifurcation point 50
N = 1; % number of bifurcation points to consider 41
tmax = 1e4; % how long to simulate (if no clear. beyond this = persistent)
repRateValue = 'default'; % do not alter
kappaStyle = 'exponential'; % default = 'exponential', 'uniform' for Fig S6D
plotOrNot = false; % whether to plot the result if want to visualize

% Fit parameters
h = load('AvC - Fits/ga-2022-10-24-at-05-13-49/fitdata.mat');

kmode_var = fliplr(logspace(2,4,N)); % varied mode values
kspan_var = linspace(1.1,3,N); % varied stdev of log-normal distribution

% Compute clearance times from varying the mode and span
data_mode = computeData(kmode_var,tmax,n,N,'kmode',[],...
    kappaStyle,repRateValue,h);
data_span = computeData(kspan_var,tmax,n,N,'stdev','same',...
    kappaStyle,repRateValue,h);

%% Plot the data

if plotOrNot
    figure(3), w = 17.8; ht = 6; % width and height, in cm
    set(gcf,'units','centimeters','position',[0 2 w ht],...
        'paperunits','centimeters','papersize',[w+eps ht+eps])

    % Plot effect of varying the mode
    subplotter(1,2,1)
    for i = 1:N
        scatter(repmat(1./kmode_var(i),[n,1])+...
            1./kmode_var(i)/200.*randn(n,1),...
            data_mode{i},4,'ko','markerfacecolor','k')
        hold on
    end
    set(gca,'xscale','log','yscale','log')
    plotter([],[],'k','','Distribution mode','Time to clear.',[],'days')
    ylim([1 1e4])
    hold off

    % Plot effect of varying the span
    subplotter(1,2,2)
    for i = 1:N
        scatter(repmat(kspan_var(i),[n,1])+...
            kspan_var/400.*randn(n,1),...
            data_span{i},4,'ko','markerfacecolor','k')
        hold on
    end
    set(gca,'yscale','log')
    plotter([],[],'k','','Distribution span','Time to clear.',[],'days')
    ylim([1 1e4])
    hold off
end

%% Function that outputs a cell for the different values per point
function [data] = computeData(...
    bfnpar,tmax,n,N,bfnparname,bounds,kappaStyle,repRateValue,h)

data = cell(N,1); % cell that will contain the clearance time data

% Load up parameters
bestpar = h.fitspecs.bestpar;
sigma_E = bestpar(1);   %
r_E = bestpar(2);       %
delta_E = bestpar(3);   %
epsilon = bestpar(4);   %
Vmax = bestpar(5);      %
if strcmp(repRateValue,'default') && strcmp(kappaStyle,'exponential')
    beta_c = bestpar(6);      % default, Fig. 3A (mode), 3B (span)
elseif strcmp(repRateValue,'reduced') && strcmp(kappaStyle,'exponential')
    beta_c = 1.16;  % not included in paper
elseif strcmp(kappaStyle,'uniform') && strcmp(repRateValue,'default')
    beta_c = 1.15; % uniform distribution, Fig. S6D (mode)
else 
    if ~strcmp(repRateValue,'default') && ~strcmp(repRateValue,'reduced')
        error(['Must set ''repRateValue'' to either ''default'' ',...
            'or ''reduced'''])
    elseif ~strcmp(kappaStyle,'uniform')&&~strcmp(kappaStyle,'exponential')
        error(['Must set ''kappaStyle'' to either ''exponential'' ',...
            'or ''uniform'''])
    else
        error(['Not a valid combination for ''repRateValue'' and ',...
            '''kappaStyle''.'])
    end
end
% beta_a = bestpar(7);    %
kappa_E = bestpar(8);   %
kappa_V = bestpar(9);   %
kmode = 10.^bestpar(10);  %
a = bestpar(11);        %
b = bestpar(12);        %
WTstdev = bestpar(13);    %
kappa_E_min = bestpar(14);
% ICa = bestpar(15)
ICc = bestpar(16);

% Loop through different mode or span values
for ii = 1:N
    estimateTimeLeft(ii,1,N,1)
    current_values = zeros(n,1);
    current_par = bfnpar(ii);

    % Select which parameter (mode or span)
    if strcmp(bfnparname,'kmode')
        kmode = current_par;
        stdev = WTstdev;
    elseif strcmp(bfnparname,'stdev')
        stdev = current_par;
    end
    
    % Use parallel computing to speed up loop for different realizations at
    % each mode/span value. If throws error, change parfor to for. It will
    % take longer.
    parfor jj = 1:n
        % Run the simulations
        [tc_tdt,solc_tdt,~] = getSimulation(tmax,sigma_E, r_E, delta_E, ...
            epsilon, Vmax, beta_c, kappa_E,kappa_V, kmode, a, b,...
            stdev,kappa_E_min,ICc,bounds,WTstdev,kappaStyle);

        % Determine the time under until virus is depleted
        tdeplete_tdt = tc_tdt(solc_tdt.V<1);
        if isempty(tdeplete_tdt)
            current_values(jj) = tmax;
        else
            current_values(jj) = tdeplete_tdt(1);
        end
    end
    data{ii} = current_values;

    clc
end

end

%% Function that generates the simulation

function [t,sol,par] = ...
    getSimulation(tmax,sigma_E, r_E, delta_E, epsilon, Vmax, beta_c, ...
    kappa_E,kappa_V, kmode, a, b, stdev, kappa_E_min,ICc,bounds,...
    WTstdev,kappaStyle)

ICs = [ICc 0.001];

if strcmp(bounds,'same')
    Amin = 10.^(-(log10(kmode)+5*log10(WTstdev)));
    Amax = 10.^(-(log10(kmode)-5*log10(WTstdev)));
else
    Amin = 10.^(-(log10(kmode)+5*log10(stdev)));
    Amax = 10.^(-(log10(kmode)-5*log10(stdev)));
end

[t, sol, par] = AvC_2_3('chronic','WT',tmax,'sigma_E',sigma_E,...
    'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Vmax,...
    'beta_c',beta_c,'kappa_E',kappa_E,'kappa_V',kappa_V,...
    'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
    'deterministic','off','ICs',ICs,'k_std',stdev,'kstyle','log',...
    'kappa_E_min',kappa_E_min,'diststyle_kappa_E',kappaStyle);

end