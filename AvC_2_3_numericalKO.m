%% This script is to determine effects of dTCR repertoires (Figs 4,S3,S6)
% Creation date: May 5th 2021, (C) Hassan (Sam) Jamaleddine

% This script takes some time to run- the results were therefore saved as
% .mat files, and loaded in AvC_figs_main or AvC_figs_supp to plot the
% results:
%
%   AvC_2_3_numericalKO.mat: time to clearance data shown in Fig. 4A-B,
%       with default replication rate obtained from parameter fitting and
%       exhaustion rate drawn from exponential distribution
%
%   AvC_2_3_numericalKO_reducedRep.mat: time to clearance data shown in
%       Fig. S3A-B, where the acute cluster can be seen at high
%       pMHC-reactivity cut-off thresholds due to more high-affinity T
%       cells
%
%   AvC_2_3_numericalKO_reducedRep_notConserved.mat: time to clearance data
%       shown in Fig. S3C-D, where the acute cluster disappears again at
%       high cut-off thresholds since there is no longer a compensatory
%       increase in high-affinity T cells
%
%   AvC_2_3_numericalKO_uniform.mat: time to clearance data shown in 
%       Fig. S6E, where exhaustion rates are instead drawn from a uniform
%       distribution.

%% Do the computation

n = 50; % number of simulations at each bifurcation point
N = 41; % number of bifurcation points to consider
tmax = 1e4; % how long to simulate (if no clear. beyond this = persistent)
getRidOf = 'low'; % do not change ('high' = remove high reactivity clones)
keepTotConstant = true; % true = Fig 4A-B and S3A/B, false = Fig S3C/D
repRateValue = 'default'; % default = Fig 4A-B, 'reduced' for all Fig S3
kappaStyle = 'exponential'; % 'exponential': Fig 4A-B/S3, 'uniform': Fig S6E
plotOrNot = false; % whether to plot results for visualization

% Load fit parameters
h = load('AvC - Fits/ga-2022-10-24-at-05-13-49/fitdata.mat');

kmode = 10.^(h.fitspecs.bestpar(10)); % pMHC-reactivity mode
kspan = h.fitspecs.bestpar(13); % pMHC-reactivity span
Amin = 10.^(-(log10(kmode)+5*log10(kspan))); % min pMHC-reactivity
Amax = 10.^(-(log10(kmode)-5*log10(kspan))); % max pMHC-reactivity
k = logspace(log10(1./Amax),log10(1./Amin),500); % vector, 1/reactivity

if strcmp(getRidOf,'high') % what happens if we get rid of high affinity
    xcutoff = logspace(log10(1./Amax),log10(kmode),N)';
    mask = heaviside(k-xcutoff);
elseif strcmp(getRidOf,'low') % what happens if we get rid of low affinity
    xcutoff = logspace(log10(kmode),log10(1./Amin),N)';
    mask = heaviside(-k+xcutoff);
end

% Compute times to pathogen clearance
[datattc, par, Enaive] = computeData(mask,tmax,n,N,kmode,...
    keepTotConstant,kappaStyle,repRateValue,h);

%% Plot results (if plotOrNot = true)

if plotOrNot
    figurer(3,'w',8.6)
    colors = repmat(logspace(log10(0.8),log10(0.2),N)',[1 3]);
    for i = 1:N
        plotter(1./k,par(i,:),colors(i,:),'','pMHC reactivity',...
            'Thymic input')
        hold on
    end
    set(gca,'XScale','log'), hold off
    xlim([min(1./k),max(1./k)])

    figurer(4,'w',8.6,'Position',[8.6 20])
    % First plot with all the data points
    % subplotter(1,2,1,'lpad',1.5)

    for i = 1:N
        scatter(repmat(1./xcutoff(i),[n,1])+...
            1./xcutoff(i)/200.*randn(n,1),...
            datattc{i},4,'ko','markerfacecolor',colors(i,:),...
            'markeredgecolor',colors(i,:))
        hold on
    end
    set(gca,'xscale','log','yscale','log')
    plotter([],[],'k','','Avidity cutoff','Time to clearance',[],'days')
    hold off
end

%% Function that outputs a cell for the different values per point
function [datattc, sigmas, Enaive] = computeData(...
    bfnpar,tmax,n,N,kmode,keepTotConstant,kappaStyle,repRateValue,h)

datattc = cell(N,1); % cell where the clearance times will be saved
sigmas = zeros(N,500); % thymic input at different cut-offs
Enaive = zeros(N,500); % initial T cell count distributions (Fig. 4A)

bestpar = h.fitspecs.bestpar; % load in the fit parameters
sigma_E = bestpar(1);   % 
r_E = bestpar(2);       % 
delta_E = bestpar(3);   % 
epsilon = bestpar(4);   % 
Vmax = bestpar(5);      % 

% Determine what replication rate to use
if strcmp(repRateValue,'default') && strcmp(kappaStyle,'exponential')
    beta_c = bestpar(6); % default everything, Fig. 4A-B
elseif strcmp(repRateValue,'reduced') && strcmp(kappaStyle,'exponential')
    beta_c = 1.16; % reduced replication rate for all Fig. S3
elseif strcmp(kappaStyle,'uniform') && strcmp(repRateValue,'default')
    beta_c = 1.15; % replication rate to use for Fig. S6E
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
if isempty(kmode)
    kmode = 10.^bestpar(10);  %
end
a = bestpar(11);        % 
b = bestpar(12);        % 
stdev = bestpar(13);    %
kappa_E_min = bestpar(14);
% ICa = bestpar(15);
ICc = bestpar(16);

% Loop over the different values of the pMHC-reactivity cut-offs
for ii = 1:N
    estimateTimeLeft(ii,1,N,1)
    current_values_ttc = zeros(n,1);
    current_par = bfnpar(ii,:);

    [~,solc_init,pc_init] = ...
        getSimulation(0.1,sigma_E, r_E, delta_E, ...
        epsilon, Vmax, beta_c, kappa_E,kappa_V, kmode, a, b,...
        stdev,kappa_E_min,ICc,current_par,keepTotConstant,kappaStyle);

    sigmas(ii,:) = pc_init.sigma_E;
    Enaive(ii,:) = solc_init.E(1,:);
    
    % Use parallel computing to speed it up. If this throws an error for
    % some reason just change 'parfor' to for. It will take longer.
    parfor jj = 1:n
        % Run the model simulations
        [tc_KO,solc_KO,~] = ...
            getSimulation(tmax,sigma_E, r_E, delta_E, ...
            epsilon, Vmax, beta_c, kappa_E,kappa_V, kmode, a, b,...
            stdev,kappa_E_min,ICc,current_par,keepTotConstant,...
            kappaStyle);
        
        % Determine the time under until pathogen is depleted
        tdeplete_tdt = tc_KO(solc_KO.V<1);
        if isempty(tdeplete_tdt)
            current_values_ttc(jj) = tmax;
        else
            current_values_ttc(jj) = tdeplete_tdt(1);
        end
    end
    datattc{ii} = current_values_ttc;    
    clc
end

end

%% Function that generates the model time series at each value

function [t,sol,par] = ...
    getSimulation(tmax,sigma_E, r_E, delta_E, epsilon, Vmax, beta_c,...
    kappa_E,kappa_V, kmode, a, b, stdev, kappa_E_min,ICc, mask,...
    conserveTot, kappaStyle)

ICs = [ICc 0.001];

Amin = 10.^(-(log10(kmode)+5*log10(stdev)));
Amax = 10.^(-(log10(kmode)-5*log10(stdev)));

[t, sol, par] = AvC_2_3('chronic','WT',tmax,'sigma_E',sigma_E,...
    'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Vmax,...
    'beta_c',beta_c,'kappa_E',kappa_E,'kappa_V',kappa_V,...
    'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
    'deterministic','off','ICs',ICs,'k_std',stdev,'kstyle','log',...
    'mask',mask,'conserveTot',conserveTot,'kappa_E_min',kappa_E_min,...
    'diststyle_kappa_E',kappaStyle);

end