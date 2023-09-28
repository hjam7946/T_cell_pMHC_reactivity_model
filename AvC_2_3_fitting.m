% This script is to generate fit predictions using the full model,
% optimized by having options for which parameters to fit etc.
%
% Date of creation: October 26th 2020.

partofit = ["sigma_E", "r_E", "delta_E", "epsilon", "Vmax", "beta"...
    "beta_a", "kappa_E", "kappa_V", "K_E", "a","b", "stdev",...
    "kappa_E_min","ICa","ICc"]; 
useParallel = true;
useCedar = false;
useBeluga = false;
forceStability = 1; % ensure stability of S0 steady state, see SI text
useLinear = 0; % only fit to log-scaled viral loads

maxGen = 200; % maximum number of generations
popSize = 250; % population size
clusterTime = 5; % alloted time on the cluster in hours

if useCedar || useBeluga % this code is suggested by Compute Canada
    % Create a "local" cluster object
    local_cluster = parcluster('local');
    
    % Modify the JobStorageLocation to $SLURM_TMPDIR
    local_cluster.JobStorageLocation = getenv('SLURM_TMPDIR');
    
    % Start the parallel pool
    parpool(local_cluster)
    fprintf("Successfully initiated parpool \n")
end

%% Load the data

load('LCMVserumWherry.mat')

time_acute = serum.acute.mod(:,1);
time_chronic = serum.chronic.mod(:,1);
virus_acute = serum.acute.mod(:,2);
virus_chronic = serum.chronic.mod(:,2);

% Initialize parameters (doesn't matter, all are sampled from range)
initpar = [31.8079 3.1793 0.2699 3.4985e-06 2.2784e+05 0.8346 0.3485...
    3.9457 0.0901 4.1176 0.0017 23.1930 4.8150 0.2 2e-4 9.13e-4];
parameters = [... % ranges indicated in cols 2 and 3 below
    initpar(1), 10, 100;            % sigma_E
    initpar(2), 0.5, 5;             % r_E
    initpar(3), 0.05, 0.5;          % delta_E
    initpar(4), 1e-7, 5e-6;         % epsilon
    initpar(5), 1e5, 5e5            % Pmax
    initpar(6), 0.7, 1.3            % rP_c
    initpar(7), 0.7, 1              % rP_a
    initpar(8), 0.1, 5              % kappa_E
    initpar(9), 0.01, 0.1           % kappa_P
    initpar(10), log10(3.5e2), log10(1e5)      % kmode, on LOG
    initpar(11), 1e-3, 0.01         % a
    initpar(12), 1,  100            % b
    initpar(13), 1.5, 5             % std
    initpar(14), 0.1, 1             % kappa_E_min
    initpar(15), 1e-4, 3e-4         % ICa (initial condition for Va)
    initpar(16), 1e-4, 2e-3];       % ICc (initial condition for Vc)

%% Fit using genetic algorithm

nvars = length(partofit);
if ~useCedar && ~useBeluga
    maxTime = 9; % max time, in hours
    maxTime = maxTime*3600; % max time, in seconds
else
    maxTime = (clusterTime*60 - 15)*60; % max time, with 15 minute pad
end


% Get upper and lower bounds
% tic
[~,lb,ub] = setPars(ones(1,nvars),partofit,parameters);
opts = optimoptions('ga',...
    'MaxGenerations',maxGen,...         % suggested 200*nvars 
    'PopulationSize',popSize,...          % suggested 200
    'PlotFcn',{@gaplotdistance,@gaplotrange},...
    'UseParallel',useParallel,...
    'MaxTime',maxTime);
[bestpar,fval,exitflag,output,population,scores] = ...
    ga(@(X)errorfunc(X,partofit,parameters,time_chronic,virus_chronic,...
    time_acute,virus_acute,forceStability,useLinear),...
    nvars,[],[],[],[],lb,ub,[],opts);

    
%% Plot

% Assign parameters for post-fit analysis
parPreds = setPars(bestpar,partofit,parameters);

sigma_E = parPreds(1); r_E = parPreds(2); delta_E = parPreds(3);
epsilon = parPreds(4); Vmax = parPreds(5); beta = parPreds(6);
beta_a = parPreds(7); kappa_E = parPreds(8); kappa_V = parPreds(9);
K_E = 10.^parPreds(10); a = parPreds(11); b = parPreds(12);
stdev = parPreds(13); kappa_E_min = parPreds(14); ICa = parPreds(15);
ICc = parPreds(16);

% Create a struct with the bestfit parameters (and the defaults if appl.)
parstruct.sigma_E = parPreds(1); parstruct.r_E = parPreds(2);
parstruct.delta_E = parPreds(3); parstruct.epsilon = parPreds(4);
parstruct.Vmax = parPreds(5); parstruct.beta = parPreds(6);
parstruct.beta_a = parPreds(7);
parstruct.kappa_E = parPreds(8); parstruct.kappa_V = parPreds(9);
parstruct.K_E = 10.^parPreds(10); parstruct.a = parPreds(11);
parstruct.b = parPreds(12); parstruct.stdev = parPreds(13);
parstruct.kappa_E_min = parPreds(14); parstruct.ICa = parPreds(15);
parstruct.ICc = parPreds(16);

Amin = 10.^(-(log10(K_E)+5*log10(stdev)));
Amax = 10.^(-(log10(K_E)-5*log10(stdev)));

% Colours
acolor = [137 059 179]./255;
ccolor = [242 155 041]./255;
acolor0 = [244 224 255]./255;
acolor1 = [235, 199, 255]./255;
acolor2 = [219, 153, 255]./255;
ccolor0 = [255 224 184]./255;
ccolor1 = [255, 225, 186]./255;
ccolor2 = [252, 200, 129]./255;

% Plotting specs
cbar_bounds = [0.7 1];
cbar_width = 0.02;
rfk = 5; % reduction factor for ak to reduce size of 3d/surf plots
rft = 1; % reduction factor for time

figurer(10,'height',10)

% Acute
subplotter(2,2,1)
[tpreda, Vpreda, pa] = getModelPrediction([0 80],parPreds,'acute');
p3 = plotter(time_acute,virus_acute+1,'ko'); hold on
set(p3,'markersize',5,'markeredgecolor','k',...
    'LineWidth',1), hold on
p4 = plotter(tpreda,Vpreda+1,acolor,...
    '','Time','Viral Load','days','pfu/mL'); hold off
set(gca,'yscale','log')

subplotter(2,2,3)
plotAffinity(tpreda,pa,rft,rfk,acolor,acolor0,acolor1,acolor2,...
    cbar_bounds,cbar_width), hold off
ylim([Amin Amax]), xlim([0 80])


% Chronic
subplotter(2,2,2)
[tpredc, Vpredc, pc] = getModelPrediction([0 80],parPreds,'chronic');
p1 = plotter(time_chronic,virus_chronic+1,'ko');
set(p1,'markersize',5,'markeredgecolor','k',...
    'LineWidth',1), hold on
p2 = plotter(tpredc,Vpredc+1,ccolor,...
    '','Time','Viral Load','days','pfu/mL'); hold off
set(gca,'yscale','log')

subplotter(2,2,4)
plotAffinity(tpredc,pc,rft,rfk,ccolor,ccolor0,ccolor1,ccolor2,...
    cbar_bounds,cbar_width), hold off
ylim([Amin Amax]), xlim([0 80])
    

%% Define the outputs (e.g. specs, best fit parameters)

% Print the parameters to the MATLAB command window
fprintf('\nBest-fit parameter results:\n')
for ii = 1:length(bestpar)
    str = "   " + partofit(ii) + ": " + num2str(bestpar(ii)) + "\n";
    fprintf(str)
end

% Create filenames for the output and output figure
if useCedar
    t1 = datetime; t1.TimeZone = '-08:00';
else
    t1 = datetime; t1.TimeZone = '-05:00';
end

t2 = datetime(t1,'TimeZone','-05:00');
theyear = char(string(year(t2)));
themonth = char(string(month(t2)));
if length(themonth) == 1
    themonth = ['0' themonth];
end
theday = char(string(day(t2)));
if length(theday) == 1
    theday = ['0' theday];
end
thehour = char(string(hour(t2)));
if length(thehour) == 1
    thehour = ['0' thehour];
end
theminute = char(string(minute(t2))); 
if length(theminute) == 1
    theminute = ['0' theminute];
end
thesecond = char(string(second(t2)));
if strcmp(thesecond(2),'.')
    thesecond = ['0' thesecond(1)];
else
    thesecond = [thesecond(1) thesecond(2)];
end
thedate = ['ga-' theyear '-' themonth '-' theday ...
    '-at-' thehour '-' theminute '-' thesecond];
if ~isfolder([cd '/ga_outputs/' thedate])
    mkdir([cd '/ga_outputs/' thedate])
    thepath = [cd '/ga_outputs/' thedate];
else
    mkdir([cd '/ga_outputs/' thedate '-v2'])
    thepath = [cd '/ga_outputs/' thedate '-v2'];
end
savename_data = [thepath '/fitdata.mat'];
savename_specfig = [thepath '/specfig.fig'];
savename_resultfig = [thepath '/resultfig.fig'];
disp(savename_resultfig)

% Create a flag for the numerical stability of the lower steady state
Ess = pc.Ess';
k = pc.k;
V0stable = sum(kappa_V.*Ess./k) > a .* beta;

% Create struct with all necessary elements
fitspecs.partofit = partofit;
fitspecs.V0stable = V0stable;
fitspecs.bestpar = bestpar;
fitspecs.fval = fval;
fitspecs.output = output;
fitspecs.exitflag = exitflag;
fitspecs.population = population;
fitspecs.scores = scores;
fitspecs.gaopts = opts;
fitspecs.useLinear = useLinear;
fitspecs.forceStability = forceStability;
fitspecs.maxGen = maxGen;
fitspecs.popSize = popSize;
fitspecs.clusterTime = clusterTime;

% Figure names
pause(30.*rand(1,1)) % decrease chance of overwriting files
fig_name = ['|| ' theyear '-' themonth '-' theday ' || '...
    thehour ':' theminute '||'];
set(figure(1),'Name',['Specs ' fig_name])
set(figure(10),'Name',['Result ' fig_name])


% Output!
save(savename_data,'fitspecs','parstruct','forceStability','useLinear')
savefig(figure(1),savename_specfig)
savefig(figure(10),savename_resultfig)

fprintf("\n")
fprintf("The MATLAB code has reached the end successfully \n")


%% Error function evaluated by the genetic algorithm
function error_tot = errorfunc(par,partofit,parameters,time_c,virus_c,...
    time_a,virus_a,forceStability,useLinear)

% par is the parameter value being evaluated at the time
% partofit is the string array containing names of parameters being fit
% parameters is the matrix containing all parameter starts, mins and maxes
% the rest is the data

partorun = setPars(par,partofit,parameters);

sigma_E = partorun(1);
r_E = partorun(2);
delta_E = partorun(3);
epsilon = partorun(4);
Vmax = partorun(5);
beta = partorun(6);
beta_a = partorun(7);
kappa_E = partorun(8);
kappa_V = partorun(9);
K_E = 10.^(partorun(10)); % defined K_E on log scale
a = partorun(11);
b = partorun(12);
stdev = partorun(13);
kappa_E_min = partorun(14);
ICa = partorun(15);
ICc = partorun(16);

% Chronic prediction
[~,Vc,pout] = AvC_2_3_getModelPrediction(time_c,sigma_E, r_E, delta_E,...
    epsilon, Vmax, beta, kappa_E, kappa_V, K_E, a, b, stdev, ...
    kappa_E_min, ICa, ICc, 'chronic');

w_c = [1 1 1 1 1 1 1 1 2 3 2]'; % weights for chronic data (ensure clear.)
if useLinear
    error_lin_c = 1e-5.^2 * sum(w_c.*(Vc - virus_c).^2)./sum(w_c); % sse
else
    error_lin_c = 0;
end
error_log_c = sum(w_c.*(log10(Vc+1) - log10(virus_c+1)).^2)./sum(w_c);
error_c = 2*error_lin_c + error_log_c;

% Acute prediction
[~,Va,~] = AvC_2_3_getModelPrediction(time_a,sigma_E, r_E, delta_E,...
    epsilon, Vmax, beta_a, kappa_E, kappa_V, K_E, a, b, stdev, ...
    kappa_E_min, ICa, ICc, 'acute');

w_a = [1 1 2 2 3 1 1 0 0]'; % weights for acute data (don't fit to 0s)
error_log_a = sum(w_a.*(log10(Va+1) - log10(virus_a+1)).^2)./sum(w_a);

% Check that the derivative of a is positive and Va is large enough
error_der_a = 0; 
if Va(2) < Va(1) || Va(3) < Va(2) || Va(3) < 1e2
    error_der_a = 25;
end
error_a = error_der_a + error_log_a;

% Total error
error_tot = error_a + error_c;

% Penalize if the lower steady state is numerically unstable
kappa_V = pout.kappa_V;
k = pout.k;
Ess = pout.Ess';
a = pout.a;
beta = pout.beta;

V0stable = sum(kappa_V.*Ess./k) > a .* beta;
if ~V0stable && forceStability
    error_tot = error_tot+25;
end

% Penalize if beta_a > beta_C
if beta < beta_a
    error_tot = error_tot+25;
end

% Penalize if std is too large using the derived condition
f = 1/100;
kmin = 100;
b = sqrt(1./(2.*log(1./f)));
a = 1./kmin.^b;
stdmax = a.*K_E.^b;
if stdev > stdmax
    error_tot = error_tot+25;
end

% Penalize if solution evolves toward bounds of allowable affinity range
k_prop = pout.k_prop;
if any(k_prop(:,end)> 0.01 .* max(k_prop,[],2))
    error_tot = error_tot+25;
end

end

%% Model prediction

function [t, V, p] = getModelPrediction(tspan,partorun,type)

sigma_E = partorun(1);
r_E = partorun(2);
delta_E = partorun(3);
epsilon = partorun(4);
Vmax = partorun(5);
beta = partorun(6);
beta_a = partorun(7);
kappa_E = partorun(8);
kappa_V = partorun(9);
K_E = 10.^(partorun(10));
a = partorun(11);
b = partorun(12);
stdev = partorun(13);
kappa_E_min = partorun(14);
ICa = partorun(15);
ICc = partorun(16);

if strcmp(type,'acute')
    beta = beta_a;
end

[t,V, p] = AvC_2_3_getModelPrediction(tspan,sigma_E, r_E, delta_E,...
    epsilon, Vmax, beta, kappa_E, kappa_V, K_E, a, b, stdev, ...
    kappa_E_min, ICa, ICc, type);

end

%% Set parameters based off which parameters are to be fit

function [partorun, lb, ub] = setPars(par,partofit,parameters)

pardefault = parameters(:,1)';
parmin = parameters(:,2)';
parmax = parameters(:,3)';

if length(par) ~= length(partofit)
    error('Inconsistent lengths for number of parameters to fit')
end

parstrings = ["sigma_E","r_E","delta_E","epsilon","Vmax","beta",...
    "beta_a","kappa_E","kappa_V","K_E","a","b","stdev","kappa_E_min",...
    "ICa","ICc"];

partorun = pardefault;
lb = zeros(size(par));
ub = zeros(size(par));

jj = 1;
for ii = 1:length(parstrings)
    if any(parstrings(ii) == partofit)
        partorun(ii) = par(jj);
        lb(jj) = parmin(ii);
        ub(jj) = parmax(ii);
        jj = jj+1;
    end
end

end

%% Plot colourmap

function plotAffinity(t,par,rft,rfk,color,color0,color1,color2,...
    cbar_bounds,cbar_width)

ak = par.ak;
k = par.k;
k_prop = par.k_prop;

[x,y]=meshgrid(t(rft:rft:end),ak(rfk:rfk:end));
% z = E(rft:rft:end,rfk:rfk:end)';
z = zeros(size(x));
for jj = 1:length(t(rft:rft:end))
    z(:,jj) = k_prop(rft.*jj,rfk:rfk:end) ...
        ./ max(k_prop(rft.*jj,rfk:rfk:end))*100;
end

% cmap_a1 = customcolormap([0 1],[acolor; 1 1 1]);
cmap = customcolormap([0 0.5 0.9 1],...
    [color2; color1; color0; [1 1 1]]);
surf(x,y,z,'edgecolor','none'); view(0,90), hold on
plot3([0 t(end)], repmat(k_prop(1,:)*(1./k),[1 2]),[109 109],...
    'Color', [0.6 0.6 0.6],'LineWidth',0.5)
ak_avg = k_prop*(1./k);
plot3(t,ak_avg,110.*ones(size(t)),'Color',color,...
    'LineWidth',1.5)
% plot3([ttc_acute ttc_acute],[Amin Amax],[110 110],'Color','k',...
%     'LineWidth',0.5,'LineStyle','--')
set(gca,'YScale','log','Layer','top',...
    'XGrid','off','YGrid','off','XMinorGrid','off',...
    'YMinorGrid','off','Colormap',cmap)
% xlim(xlimacute)
plotter([],[],'k','','\rmTime','\rmpMHC reactivity','days',[],...
    'AxFontSize',7)
ax = gca; ax_pos = get(ax,'Position');
cbar = colorbar(ax,'TickDirection','out','Box','off',...
    'Limits',[0 100],'Ticks',[0 100],'Location',...
    'northoutside','AxisLocation','in');
ax.Position = ax_pos;
cbar_pos = cbar.Position;
cbar_xpos = [cbar_pos(1) cbar_pos(1)+cbar_pos(3)];
cbar_length = diff(cbar_xpos).*diff(cbar_bounds);
cbar_lbound = cbar_pos(1)+cbar_bounds(1).*diff(cbar_xpos);
cbar.Position = [cbar_lbound cbar_pos(2) cbar_length cbar_width];
cbar.TickLength = 0.07;

end