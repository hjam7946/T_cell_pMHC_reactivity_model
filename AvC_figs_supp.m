% Run this code to generate the modelling parameters and results shown in
% the figures of the main manuscript (Figures S1A-C, S2, S3, S5 and S6) as
% well as the supplemental movies (Movies S1 and S2)
%
% Note that, since some results are generated from time-consuming
% simulations, these are instead loaded from .mat files. See where and how
% to regenerate them, on a case-by-case basis, below.
%
% (C) Hassan (Sam) Jamaleddine, 2022

%% Options

fig_num = 'S6'; % choose what to plot/run (S1,S2,S3,S5,M1,M2)

% Load parameters from fit result and simulate
load('AvC - Fits/ga-2022-10-24-at-05-13-49/fitdata.mat');
bestpar=fitspecs.bestpar; % Parameter outputs from genetic algorithm

sigma_E = bestpar(1);       % total thymic input
r_E = bestpar(2);           % T-cell activation/proliferation rate
delta_E = bestpar(3);       % natural t-cell turnover
epsilon = bestpar(4);       % homeostatic deletion rate
Pmax = bestpar(5);          % maximum viral load (Pmax in paper)
rP_c = bestpar(6);          % chronic viral replication rate
rP_a = bestpar(7);          % acute viral replication rate
kappa_E = bestpar(8);       % avg pathogen-dependent T cell exhaustion
kappa_P = bestpar(9);       % pathogen clearance rate by T cells
kmode = 10.^bestpar(10);    % mode value of 1/reactivity (mean of logs)
a = bestpar(11);            % scaling factor in pathogen clearance
b = bestpar(12);            % scaling factor in t cell depletion
kspan = bestpar(13);        % span of the pMHC reactivity distribution
kappa_E_min = bestpar(14);  % minimum exhaustion rate
ICa = bestpar(15);          % initial pathogen load for acute
ICc = bestpar(16);          % initial pathogen load for chronic

Amin = 10.^(-(log10(kmode)+5*log10(kspan))); % define min reactivity
Amax = 10.^(-(log10(kmode)-5*log10(kspan))); % define max reactivity
kstyle = 'log'; % log-normally distributed reactivity (vs. gaussian)
direction = 'descend'; % do not change

%% Figure S1 - Supplemental to Fig 2

if strcmp(fig_num,'S1')
    w = 11.4;   % width in cm
    h = 4.5;    % height in cm
    ms = 5;    % marker size, in pts
    hratios = [0.38 0.38 0.24];
    computettcs = false; % set to true to recompute clearance times
    tstep = 0.1;
    sim_t_acute = (0:tstep:100)';
    sim_t_chronic = (0:tstep:100)';

    close(figure(2))
    figurer(2,'w',w,'h',h,'PaperSize','letter')

    % Colors
    acolor = [137 059 179]./255;
    ccolor = [242 155 041]./255;

    acolor0 = [244 224 255]./255;
    acolor1 = [235, 199, 255]./255;
    acolor2 = [219, 153, 255]./255;
    ccolor0 = [255 224 184]./255;
    ccolor1 = [255, 225, 186]./255;
    ccolor2 = [252, 200, 129]./255;

    % Load and process the serum data from Wherry et al., 2003
    load('LCMVserumWherry.mat')
    data_t_acute = serum.acute.mod(:,1);
    data_V_acute = serum.acute.mod(:,2);
    data_t_chronic = serum.chronic.mod(:,1);
    data_V_chronic = serum.chronic.mod(:,2);

    deterministic = 'on'; % remove effect of pramaeter randomization

    % Uncomment if want to run simulation

    [~, sol_acute, ~] = ...
        AvC_2_3('acute','WT',sim_t_acute,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'deterministic',deterministic,...
        'ICs',[ICa,0.001],'kappa_E_min',kappa_E_min);
    [~, sol_chronic,~] = ...
        AvC_2_3('chronic','WT',sim_t_chronic,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',2,'kstyle',kstyle,'deterministic',...
        deterministic, 'ICs', [ICc 0.001], 'kappa_E_min',kappa_E_min);

    sim_V_acute = sol_acute.V;
    sim_V_chronic = sol_chronic.V;

    % Acute (data + simulation)
    subplotter(1,3,1,'hratios',hratios)
    p1 = plotter(sim_t_acute,sim_V_acute+1,acolor); hold on
    p2 = plotter(data_t_acute,data_V_acute+1,'o','A','\rmTime',...
        '\rmViral load + 1','days','PFU/mL','AxFontSize',7); hold off

    set(p2,'MarkerSize',ms,'MarkerEdgeColor',...
        'k','MarkerFaceColor','none','LineWidth',0.5,'Color','none')
    set(gca,'YScale','log','YLim',[1 1e6],'YTick',[1 1e2 1e4 1e6],...
        'YMinorTick','off')

    % Chronic (data + simulation)
    subplotter(1,3,2,'hratios',hratios)
    p3 = plotter(sim_t_chronic,sim_V_chronic+1,ccolor); hold on
    p4 = plotter(data_t_chronic,data_V_chronic+1,'o','B','\rmTime',...
        '\rmViral load + 1','days','PFU/mL','AxFontSize',7); hold off

    set(p4,'MarkerSize',ms,'MarkerEdgeColor',...
        'k','MarkerFaceColor','none','LineWidth',0.5,'Color','none')

    set(gca,'YScale','log','YLim',[1 1e6],'YTick',[1 1e2 1e4 1e6],...
        'YMinorTick','off')
    xlim([0 100])

    % PANEL C - (Compute and) plot times to clearance for AvC

    % Define and compute times to clearance for acute and chronic
    n = 100; % number of simulations per condition (acute or chronic)
    tmax = 1e3;
    
    % Compute or load clearance times
    if computettcs
        ttc_acute = zeros(1,n);
        ttc_chronic = zeros(1,n);
        parfor ii = 1:n
            % estimateTimeLeft(ii,1,n,1)
            [ta,sola,~] = ... % simulation for acute
                AvC_2_3('acute','WT',200,'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_a',rP_a,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'kstyle',kstyle,...
                'deterministic','off','ICs',[ICa,0.001],...
                'kappa_E_min',kappa_E_min);
            [tc,solc,~] = ... % simulation for chronic
                AvC_2_3('chronic','WT',200,'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_c',rP_c,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'kstyle',kstyle,...
                'deterministic','off','ICs',[ICc,0.001],...
                'kappa_E_min',kappa_E_min);

            % Compute clearance times
            ttc_acute(ii) = computeTTC(ta,sola.V);
            ttc_chronic(ii) = computeTTC(tc,solc.V,tmax);

            fprintf(num2str(ii) + "\n")
        end
    else
        ttcdata = load('AvC_2_3_defaultTTCs.mat');
        ttc_acute = ttcdata.ttc_acute;
        ttc_chronic = ttcdata.ttc_chronic;
    end

    % Plot the times to clearance for acute or for chronic as bar plot
    subplotter(1,3,3,'hratios',hratios)
    bwidth = 0.6;

    ttc_acute = sort(ttc_acute);
    ttc_chronic = sort(ttc_chronic);

    % Bar graph
    br = bar([1,2],[median(ttc_acute);median(ttc_chronic)],bwidth);
    br.FaceColor = 'flat';
    br.CData(1,:)=acolor2;
    br.CData(2,:)=ccolor2;
    br.LineWidth = 0.75; hold on

    mu_a = median(ttc_acute);
    mu_c = median(ttc_chronic);

    neg_a = mu_a - min(ttc_acute(5:95));
    pos_a = max(ttc_acute(5:95))-mu_a;
    neg_c = mu_c - min(ttc_chronic(5:95));
    pos_c = max(ttc_chronic(5:95))-mu_c;

    % Error bars
    hold on
    errorbar([1,2],[mu_a mu_c],[neg_a neg_c],[pos_a pos_c],...
        'marker','none','linestyle','none','Color','k',...
        'Capsize',5,'LineWidth',1)
    set(gca,'YScale','log','XTickLabels',{'Arm','Cl13'})
    plotter([],[],'k','C','','\rmTime to clearance',[],'days')
    ylim([1 100]), xlim([0.2 2.8])
    hold off

    [~,pvalue] = ttest2(ttc_acute,ttc_chronic,'Vartype','unequal');
    sigstar([1,2],pvalue);
end

%% Figure S2 - Plot of the bifurcation analysis of the system

if strcmp(fig_num,'S2')

    w = 17.4;   % width in cm
    h = 5;    % height in cm
    kappa_lim = [0 4];
    ak_lim = [1e-5 1e-1];

    close(figure(3))
    figurer(3,'w',w,'h',h,'PaperSize','letter')
    
    % Load the bifurcation data from XPP_AUTO
    addpath([cd '\AvC - Bifurcation Data\'])
    
    % Load and process the 1par V vs kappa_E bfn data
    datakexp = load('AvC_2_3_logVvsKexp_revised.dat');
    datakexp(:,1) = 10.^(-datakexp(:,1));
    datakexp(:,[2 3]) = 10.^(datakexp(:,[2 3]));
    datakexp = datakexp(10.^(datakexp(:,2))>0,:);

    % Load and process the 1par V vs kappa_E bfn data
    datakappa_E = load('AvC_2_3_logVvsKappaE_revised.dat');
    datakappa_E(:,[2 3]) = 10.^(datakappa_E(:,[2 3]));
    datakappa_E = datakappa_E(10.^(datakappa_E(:,2))>0,:);

    % Load and process the 2par bifurcation data
    data2par = load('AvC_2_3_KexpVsKappaE_revised.dat');
    av_2par = 1./(10.^(data2par(:,2)));
    kappaE_2par = data2par(:,1);

    deterministic = 'on';

    % Run simulations for acute and chronic, WT
    [t,sol,par] = ...
        AvC_2_3('chronic','WT',200,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',2,'kstyle',kstyle,'deterministic',...
        deterministic,'kappa_E_min',kappa_E_min);
    kappE_eff = sol.k_prop*par.kappa_E;

    % Plot the 1 parameter bifurcation w.r.t. avidity
    subplotter(1,3,1) % 1 par plot
    xpp_plot('A','\rmpMHC reactivity','\rmPathogen load','\it{a_k}\rm',...
        '\it{P}\rm{ + 1}','data',datakexp)
    set(gca,'XScale','log','YScale','log','XTick',[1e-5 1e-3 1e-1])
    xlim([1e-5-eps 1e-1]), ylim([1e-1 1e6])

    % Plot the 1 parameter bifurcation w.r.t. kappa_E
    subplotter(1,3,2) % 1 par plot
    xpp_plot('B','\rmExhaustion rate','\rmPathogen load',...
        ['\it{' char(954) '_E}\rm'],...
        '\it{P}\rm{ + 1}','data',datakappa_E)
    set(gca,'YScale','log')
    xlim(kappa_lim), ylim([1e-1 1e6])

    % Plot the 2 parameter bifurcation with shaded region and trajectory
    subplotter(1,3,3) % 2 par plot
    ymax = kappa_lim(end);
    xlim([-5 -1]), ylim([0 ymax])
    plotcurvearrows(log10(1./sol.k_avg),kappE_eff,6,...
        'HeadSize',4,'HeadStyle','cback3') % arrows on curve
    fill([av_2par;min(av_2par);1e-5;1e-5;min(av_2par);0.1;0.1],...
        [kappaE_2par;0;0;ymax;ymax;ymax;kappaE_2par(1)],...
        [0.9 0.9 0.9],'edgecolor','none'), hold on % filled area
    plotter(av_2par,kappaE_2par,[0.5 0.5 0.5]) % bifurcation curve
    scatter(1./sol.k_avg(1),kappE_eff(1),13,...
        'markeredgecolor','k','markerfacecolor','k') % start point
    p=plotter(1./sol.k_avg,kappE_eff,'k','C','\rmpMHC reactivity',...
        '\rmExhaustion rate','\it{a_k}\rm',['\it{' char(954) '_E}\rm']);
    hold off
    set(gca,'XScale','log'); p.LineWidth = 1; % trajectory
    xlim([1e-5-eps 1e-1]), ylim(kappa_lim)
    set(gca,'XTick',[1e-5 1e-3 1e-1],'YTick',0:kappa_lim(end))
end

%% Figure S3 - Numerical KO w/o compensatory increase of high-avs

if strcmp(fig_num,'S3')

    w = 17.4;   % width in cm
    h = 10;    % height in cm
    tdt_col = [171, 36, 41]./255; % matcol(7,:);
    wt_col = [0.6 0.6 0.6];

    close(figure(4))
    figurer(4,'w',w,'h',h,'PaperSize','letter')

    % Numerical KO, reduced replication rate, total conserved
    % (to generate data, run AvC_2_3_numericalKO.m with repRateValue set to
    % 'reduced')
    numKO1 = load('AvC_2_3_numericalKO_reducedRep.mat');
    xcutoff1 = numKO1.xcutoff;
    datattc1 = numKO1.datattc;
    N_KO1 = numKO1.N;
    n_KO1 = numKO1.n;
    k = numKO1.k;
    Enaive1 = numKO1.Enaive;

    % ------------------------------------------------------------------ %
    % B - Numerical KO, data
    subplotter(2,2,2)
        colors = zeros(N_KO1,3);
    for ii = 1:N_KO1
        colors(ii,:) = round((wt_col-tdt_col)*255*(ii-1)./N_KO1 ...
            + tdt_col.*255)./255;
        staggerscatter([nan(1,ii-1) 1./xcutoff1(ii) nan(1,N_KO1-ii)],...
            [nan(n_KO1,ii-1) datattc1{ii} nan(n_KO1,N_KO1-ii)],4,...
            'markeredgecolor',colors(ii,:),...
            'markerfacecolor',colors(ii,:),'XScale','log',...
            'YScale','log','XSpread',0.5), hold on
        hold on
    end
    set(gca,'xscale','log','yscale','log','ytick',[1e1 1e2 1e3 1e4])
    plotter([],[],'k','B','\rmpMHC reactivity cut-off threshold',...
        '\rmTime to clear.',[],'days','AxFontSize',7)
    fs = get(get(gca,'XLabel'),'FontSize'); % carry forward for 3D plot
    ylim([1 1e4])
    hold off

    % ------------------------------------------------------------------ %
    % A - Avidity cutoff profile
    subplotter(2,2,1)
    for i = 1:N_KO1
        p = plot(1./k, Enaive1(i,:));
        p.Color = colors(i,:);
        p.LineWidth = 1;
        hold on
        if i == 1
            plotter([],[],'k','A')
        end
    end
    set(gca,'XScale','log','YScale','linear'), hold off
    xlim([min(1./k),max(1./k)])
    set(gca,'YTick',[0 1 2])
    xlabel('\rm{pMHC reactivity}','FontSize',fs)
    ylabel('\rmT cell count','FontSize',fs)
    hold off

    % Numerical KO, reduced replication rate, total NOT conserved
    % (to generate data, run AvC_2_3_numericalKO.m with repRateValue set to
    % 'reduced' and totConserved set to false)
    numKO2 = load('AvC_2_3_numericalKO_reducedRep_notConserved.mat');
    xcutoff2 = numKO2.xcutoff;
    datattc2 = numKO2.datattc;
    N_KO2 = numKO2.N;
    n_KO2 = numKO2.n;
    Enaive2 = numKO2.Enaive;

    % ------------------------------------------------------------------ %
    % B - Numerical KO, data
    subplotter(2,2,4)
        colors = zeros(N_KO2,3);
    for ii = 1:N_KO2
        colors(ii,:) = round((wt_col-tdt_col)*255*(ii-1)./N_KO2 ...
            + tdt_col.*255)./255;
        staggerscatter([nan(1,ii-1) 1./xcutoff2(ii) nan(1,N_KO2-ii)],...
            [nan(n_KO2,ii-1) datattc2{ii} nan(n_KO2,N_KO2-ii)],4,...
            'markeredgecolor',colors(ii,:),...
            'markerfacecolor',colors(ii,:),'XScale','log',...
            'YScale','log','XSpread',0.5), hold on
        hold on
    end
    set(gca,'xscale','log','yscale','log','ytick',[1e1 1e2 1e3 1e4])
    plotter([],[],'k','D','\rmpMHC reactivity cut-off threshold',...
        '\rmTime to clear.',[],'days','AxFontSize',7)
    fs = get(get(gca,'XLabel'),'FontSize'); % carry forward for 3D plot
    ylim([1 1e4])
    hold off

    % ------------------------------------------------------------------ %
    % A - Avidity cutoff profile
    subplotter(2,2,3)
    for i = 1:N_KO2
        p = plot(1./k, Enaive2(i,:));
        p.Color = colors(i,:);
        p.LineWidth = 1;
        hold on
        if i == 1
            plotter([],[],'k','C')
        end
    end
    set(gca,'XScale','log','YScale','linear'), hold off
    xlim([min(1./k),max(1./k)])
    set(gca,'YTick',[0 1 2])
    xlabel('\rm{pMHC reactivity}','FontSize',fs)
    ylabel('\rmT cell count','FontSize',fs)
    hold off

end

%% Figure S5 - Other tests of TdT roles

if strcmp(fig_num,'S5')
    w = 17.4;   % width in cm
    h = 10;    % height in cm changed from 9
    tdt_col = [171, 36, 41]./255; % matcol(7,:);
    wt_col = [0.6 0.6 0.6];
    tmax_acute = 30;
    tmax_chronic = 100;
    ms = 3; % marker size, in pts, for scatter plot
    computettcs = false; % change to 'true' to recompute from scratch
    hratio = [0.45 0.275 0.275];

    close(figure(8))
    figurer(8,'w',w,'h',h,'PaperSize','letter')
    load('colours.mat')

    % Set parameters and options for running simulations
    deterministic = 'on';
    sigma_E_tdt = sigma_E/10; % smaller precursor size in TdT KO
    n_tdt = 50; % smaller number of affinity clones
    ICa_default = ICa; % default I.C. (prior to adding random noise)
    ICc_default = ICc; % default I.C. (prior to adding random noise)

    % Run multiple times to get a distribution of times to clearance
    if computettcs == 1
        Nsims = 50; % number of simulations

        % Initialize clearance times
        ttc_wt_acute = zeros(Nsims,1);
        ttc_dTCR1_acute = zeros(Nsims,1);
        ttc_dTCR2_acute = zeros(Nsims,1);
        ttc_wt_chronic = zeros(Nsims,1);
        ttc_dTCR1_chronic = zeros(Nsims,1);
        ttc_dTCR2_chronic = zeros(Nsims,1);

        % Add variability specifically to ICs
        ICr = 0.1; % +/- 15% in log space
        ICa = 10.^(log10(ICa).*(1+ICr*randn(Nsims,1)));
        ICc = 10.^(log10(ICc).*(1+ICr*randn(Nsims,1)));

        parfor ii = 1:Nsims
            % estimateTimeLeft(ii,1,Nsims,1)
            tmax = 1e3;

            % Acute
            [ta_wt_N, sola_wt_N, ~] = ...
                AvC_2_3('acute','WT',tmax_acute,'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_a',rP_a,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',2,...
                'kstyle',kstyle,'deterministic','off',...
                'direction',direction,'ICs',[ICa(ii) 0.001],...
                'kappa_E_min',kappa_E_min);
            [ta_dTCR1_N, sola_dTCR1_N, pa_dTCR1_N] = ...
                AvC_2_3('acute','WT',tmax_acute,'sigma_E',sigma_E_tdt,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_a',rP_a,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',2,...
                'kstyle',kstyle,'deterministic','off',...
                'direction',direction,'ICs',...
                [ICa(ii) 0.001],'kappa_E_min',kappa_E_min);

            % Chronic
            [tc_wt_N, solc_wt_N,~] = ...
                AvC_2_3('chronic','WT',tmax,'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_c',rP_c,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',2,...
                'kstyle',kstyle,'deterministic','off',...
                'direction',direction,'ICs',[ICc(ii) 0.001],...
                'kappa_E_min',kappa_E_min);
            [tc_dTCR1_N, solc_dTCR1_N,pc_dTCR1_N] = ...
                AvC_2_3('chronic','WT',tmax,'sigma_E',sigma_E_tdt,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_c',rP_c,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',2,...
                'kstyle',kstyle,'deterministic','off',...
                'direction',direction,'ICs',...
                [ICc(ii) 0.001],'kappa_E_min',kappa_E_min);

            % dTCR2 (reduced N)
            % need to interpolate independent realizations of kappa_E
            % values, given how it is sampled, to not alter robustness to
            % exhaustion between the WT and the dTCR models
            x_old = linspace(log10(Amin),log10(Amax),500)';
            x_new = linspace(log10(Amin),log10(Amax),n_tdt)';
            kappaE_dTCR_a = pc_dTCR1_N.kappa_E'; % interp. INDEPENDENTLY
            kappaE_dTCR_a = interp1(x_old,kappaE_dTCR_a,x_new,'nearest');
            kappaE_dTCR_c = pa_dTCR1_N.kappa_E';
            kappaE_dTCR_c = interp1(x_old,kappaE_dTCR_c,x_new,'nearest');

            [ta_dTCR2_N, sola_dTCR2_N, ~] = ...
                AvC_2_3('acute','WT',tmax,'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_a',rP_a,'kappa_E',kappaE_dTCR_a,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',2,...
                'kstyle',kstyle,'deterministic','off',...
                'direction',direction,'ICs',[ICa(ii) 0.001],...
                'n',n_tdt,'kappa_E_min',kappa_E_min,...
                'isdist_kappa_E','off');

            [tc_dTCR2_N, solc_dTCR2_N,~] = ...
                AvC_2_3('chronic','WT',tmax,'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_c',rP_c,'kappa_E',kappaE_dTCR_c,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',2,...
                'kstyle',kstyle,'deterministic','off',...
                'direction',direction,'ICs',[ICc(ii) 0.001],...
                'n',n_tdt,'kappa_E_min',kappa_E_min,...
                'isdist_kappa_E','off');            

            % Times to clearance
            ttc_wt_acute(ii) = ...
                computeTTC(ta_wt_N,sola_wt_N.V,tmax_acute);
            ttc_dTCR1_acute(ii) = ...
                computeTTC(ta_dTCR1_N,sola_dTCR1_N.V,tmax_acute);
            ttc_dTCR2_acute(ii) = ...
                computeTTC(ta_dTCR2_N,sola_dTCR2_N.V,tmax_acute);
            ttc_wt_chronic(ii) = ...
                computeTTC(tc_wt_N,solc_wt_N.V,tmax_chronic);
            ttc_dTCR1_chronic(ii) = ...
                computeTTC(tc_dTCR1_N,solc_dTCR1_N.V,tmax_chronic);
            ttc_dTCR2_chronic(ii) = ...
                computeTTC(tc_dTCR2_N,solc_dTCR2_N.V,tmax_chronic);
        end
    else
        % Load the data that was computed previously
        % (to recompute, set computettcs = true)
        WTvsdTCR= load('AvC_2_3_WTvsdTCR_lcmv_alternative.mat');
        ttc_wt_chronic = WTvsdTCR.ttc_wt_chronic;
        ttc_dTCR1_chronic = WTvsdTCR.ttc_dTCR1_chronic;
        ttc_dTCR2_chronic = WTvsdTCR.ttc_dTCR2_chronic;
        ttc_wt_acute = WTvsdTCR.ttc_wt_acute;
        ttc_dTCR1_acute = WTvsdTCR.ttc_dTCR1_acute;
        ttc_dTCR2_acute = WTvsdTCR.ttc_dTCR2_acute;
        Nsims = length(ttc_wt_chronic);
    end

    % Run simulations for acute and chronic, WT
    [ta_wt, sola_wt,pa_wt] = ...
        AvC_2_3('acute','WT',tmax_acute,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'ICs',[ICa_default,0.001],...
        'deterministic',deterministic,'direction',direction,...
        'kappa_E_min',kappa_E_min);
    [tc_wt, solc_wt,pc_wt] = ...
        AvC_2_3('chronic','WT',tmax_chronic,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',2,'kstyle',kstyle,'deterministic',...
        deterministic,'direction',direction,'ICs',[ICc_default 0.001],...
        'kappa_E_min',kappa_E_min);

    % Run simulations for acute and chronic, TdT (reduced sigma)
    [ta_tdt1, sola_tdt1,pa_tdt1] = ...
        AvC_2_3('acute','WT',tmax_acute,'sigma_E',sigma_E_tdt,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'ICs',[ICa_default,0.001],...
        'deterministic',deterministic,'direction',direction,...
        'kappa_E_min',kappa_E_min);
    [tc_tdt1, solc_tdt1,pc_tdt1] = ...
        AvC_2_3('chronic','WT',tmax_chronic,'sigma_E',sigma_E_tdt,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',2,'kstyle',kstyle,'deterministic',...
        deterministic,'direction',direction,'ICs',[ICc_default 0.001],...
        'kappa_E_min',kappa_E_min);

    % y_out = interp1(t_out,y_out,t_interp,'spline')
    % Run simulations for acute and chronic, TdT (reduced N)
    % First, interpolate randomly sampled exhaustion rates so they are
    % distributed the same
    x_old = linspace(log10(Amin),log10(Amax),500)';
    x_new = linspace(log10(Amin),log10(Amax),n_tdt)';
    kappaE_dTCR_a = pc_tdt1.kappa_E';
    kappaE_dTCR_a = interp1(x_old,kappaE_dTCR_a,x_new,'nearest');
    kappaE_dTCR_c = pa_tdt1.kappa_E';
    kappaE_dTCR_c = interp1(x_old,kappaE_dTCR_c,x_new,'nearest');

    [ta_tdt2, sola_tdt2,pa_tdt2] = ...
        AvC_2_3('acute','WT',tmax_acute,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappaE_dTCR_a,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'ICs',[ICa_default 0.001],...
        'deterministic',deterministic,'direction',direction,...
        'n',n_tdt,'isdist_kappa_E','off','kappa_E_min',kappa_E_min);
    [tc_tdt2, solc_tdt2,pc_tdt2] = ...
        AvC_2_3('chronic','WT',tmax_chronic,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappaE_dTCR_c,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',2,'kstyle',kstyle,'deterministic',...
        deterministic,'direction',direction,'n',n_tdt,...
        'isdist_kappa_E','off','kappa_E_min',kappa_E_min,...
        'ICs', [ICc_default 0.001]);
    
    lpad = 1;
    rpad = 0.5;
    dpad = 1;
    upad = 0.65;

    % ------------------------------------------------------------------ %
    % A - Distribution, simulation (reduced N)
    subplotter(2,3,1,'hratios',hratio)
    Enaive_wt = solc_wt.E(1,:)./max(solc_wt.E(1,:))*100;
    Enaive_tdt2 = solc_tdt2.E(1,:)./max(solc_tdt2.E(1,:))*100;
    p1 = plotter(pa_wt.ak,Enaive_wt,wt_col,'A',...
        '\rmpMHC reactivity','\rmT cell count',[],'% of max'); hold on
    p2 = stem(pa_tdt2.ak,Enaive_tdt2,'Color',tdt_col,...
        'MarkerSize',ms,'MarkerFaceColor',tdt_col); hold off
    p2.BaseLine.LineStyle = 'none'; 
    set(gca,'xscale','log','XLim',[Amin Amax],...
        'XTick',[1e-5 1e-4 1e-3 1e-2],'FontSize',7)
    l = legend([p1, p2], {'WT', [char(916) 'TCR rep.']},'box','off');
    l.ItemTokenSize = [10 18];

    % ------------------------------------------------------------------ %
    % B,top - Simulation, acute
    subplotter(4,3,2,'hratios',hratio,'dpad',dpad-0.5)
    p1 = plotter(ta_tdt2,sola_tdt2.V,tdt_col); hold on
    p2 = plotter(ta_wt,sola_wt.V,wt_col,'B',''); hold off
    set(get(gca,'YLabel'),'String','Pathogen load (PFU/mL)')
    l = legend([p2 p1],{'WT',[char(916) 'TCR rep.']},'box','off');
    l.ItemTokenSize = [10 18];
    set(gca,'FontSize',7,'YScale','log','YLim',[1 1e4],...
        'YMinorTick','off')
    xlim([0 tmax_acute])

    % ------------------------------------------------------------------ %
    % B,bottom - Simulation, chronic
    subplotter(4,3,5,'hratios',hratio,'upad',upad-0.5)
    plotter(tc_tdt2,solc_tdt2.V,tdt_col); hold on
    plotter(tc_wt,solc_wt.V,wt_col,'','\rmTime','','days',[]); hold off
    xlim([0 tmax_chronic])
    set(gca,'FontSize',7,'YScale','log','YLim',[1 1e6],...
        'YMinorTick','off','YTick',[1 1e2 1e4 1e6])

    % ------------------------------------------------------------------ %
    % C - Time to clearance, WT vs. TdT KO
    
    % Acute violin plot
    subplotter(2,4,3,'hratios',hratio2,'rpad',0.4-0.4) %(2,3,6)
    vp = violinplot([ttc_wt_acute; ttc_dTCR2_acute],...
        [repmat({'WT1'},[Nsims 1]);...
        repmat({[char(916) 'TCR rep2.']},[Nsims 1])],'GroupOrder',...
        {'WT1',[char(916) 'TCR rep2.']},'ShowData',...
        true);
    vp(1).EdgeColor = wt_col; vp(2).EdgeColor = tdt_col;
    vp(1).ViolinColor = wt_col; vp(2).ViolinColor = tdt_col;
    for ii = 1:2
        vp(ii).MedianPlot.Visible = 'off';
        vp(ii).MeanPlot.Visible = 'on';
        vp(ii).BoxPlot.Visible = 'off';
        vp(ii).WhiskerPlot.Visible = 'off';
        vp(ii).ViolinPlot.LineWidth = 1;
        vp(ii).ScatterPlot.Visible = 'off';
    end
    hold off

    plotter([],[],'k','C','','\rmTime to clearance',[],'days',...
        'AxFontSize',7)
    set(get(gca,'XAxis'),'TickLabelRotation',0)
    set(gca,'XTickLabels',{'WT', [char(916) 'TCR']})
    set(gca,'YScale','linear')
    xlim([0.4 2.6]), ylim([0 20])
    
    % Chronic violin plot
    subplotter(2,4,4,'hratios',hratio2,'lpad',1.1-0.4) %(2,3,6)
    vp = violinplot([ttc_wt_chronic; ttc_dTCR2_chronic],...
        [repmat({'WT1'},[Nsims 1]);...
        repmat({[char(916) 'TCR rep1.']},[Nsims 1])],'GroupOrder',...
        {'WT1',[char(916) 'TCR rep1.']},'ShowData',...
        true);
    vp(1).EdgeColor = wt_col; vp(2).EdgeColor = tdt_col;
    vp(1).ViolinColor = wt_col; vp(2).ViolinColor = tdt_col;
    for ii = 1:2
        vp(ii).MedianPlot.Visible = 'off';
        vp(ii).MeanPlot.Visible = 'on';
        vp(ii).BoxPlot.Visible = 'off';
        vp(ii).WhiskerPlot.Visible = 'off';
        vp(ii).ViolinPlot.LineWidth = 1;
        vp(ii).ScatterPlot.Visible = 'off';
    end
    hold off

    plotter([],[],'k','','','',[],[],'AxFontSize',7)
    set(get(gca,'XAxis'),'TickLabelRotation',0)
    set(gca,'XTickLabels',{'WT', [char(916) 'TCR']})
    set(gca,'YScale','linear','YTick',[0 25 50 75 100 125])
    xlim([0.4 2.6]), ylim([0 125]),
    
    % Compute p-value + plot significance stars
    pvalue_acute = ranksum(ttc_wt_acute,ttc_dTCR2_acute);
    pvalue_chronic = ranksum(ttc_wt_chronic,ttc_dTCR2_chronic);
    disp("Rank sum p-value (acute):   " + num2str(pvalue_acute))
    disp("Rank sum p-value (chronic): " + num2str(pvalue_chronic))

    % ------------------------------------------------------------------ %
    % D - Distribution, simulation (reduced sigma)v
    subplotter(2,3,4,'hratios',hratio)
    % subplotter(2,3,3)
    Enaive_wt = solc_wt.E(1,:)./max(solc_wt.E(1,:))*100;
    Enaive_tdt1 = solc_tdt1.E(1,:)./max(solc_wt.E(1,:))*100;
    p2 = plotter(pa_tdt1.ak,Enaive_tdt1,tdt_col); hold on
    p1 = plotter(pa_wt.ak,Enaive_wt,wt_col,'D',...
        '\rmpMHC reactivity','\rmT cell count',[],'% of WT max'); hold off
    set(gca,'xscale','log','XLim',[Amin Amax],...
        'XTick',[1e-5 1e-4 1e-3 1e-2],'FontSize',7)
    l = legend([p1, p2], {'WT', [char(916) 'TCR rep.']},'box','off');
    l.ItemTokenSize = [10 18];

    % ------------------------------------------------------------------ %
    % E,top - Simulation, acute
    subplotter(4,3,8,'hratios',hratio,'dpad',dpad-0.5)
    p1 = plotter(ta_tdt1,sola_tdt1.V,tdt_col); hold on
    p2 = plotter(ta_wt,sola_wt.V,wt_col,'E',''); hold off
    set(get(gca,'YLabel'),'String','Pathogen load (PFU/mL)')
    l = legend([p2 p1],{'WT',[char(916) 'TCR rep.']},'box','off');
    l.ItemTokenSize = [10 18];
    set(gca,'FontSize',7,'YScale','log','YLim',[1 1e4],...
        'YMinorTick','off')
    xlim([0 tmax_acute])

    % ------------------------------------------------------------------ %
    % E,bottom - Simulation, chronic
    subplotter(4,3,11,'hratios',hratio,'upad',upad-0.5)
    plotter(tc_tdt1,solc_tdt1.V,tdt_col); hold on
    plotter(tc_wt,solc_wt.V,wt_col,'','\rmTime',...
        '','days',[]);
    hold off
    xlim([0 tmax_chronic])
    set(gca,'FontSize',7,'YScale','log','YLim',[1 1e6],...
        'YMinorTick','off','YTick',[1 1e2 1e4 1e6])

    % ------------------------------------------------------------------ %
    % F - Time to clearance, WT vs. TdT KO
    hratio2 = [hratio 0];
    hratio2(end-1:end) = hratio2(end-1)/2;
    Nsims = length(ttc_wt_chronic);
    addpath Violinplot-Matlab-master/
    
    % Acute violin plot
    subplotter(2,4,7,'hratios',hratio2,'rpad',0.4-0.4) %(2,3,6)
    vp = violinplot([ttc_wt_acute; ttc_dTCR1_acute],...
        [repmat({'WT1'},[Nsims 1]);...
        repmat({[char(916) 'TCR rep1.']},[Nsims 1])],'GroupOrder',...
        {'WT1',[char(916) 'TCR rep1.']},'ShowData',...
        true);
    vp(1).EdgeColor = wt_col; vp(2).EdgeColor = tdt_col;
    vp(1).ViolinColor = wt_col; vp(2).ViolinColor = tdt_col;
    for ii = 1:2
        vp(ii).MedianPlot.Visible = 'off';
        vp(ii).MeanPlot.Visible = 'on';
        vp(ii).BoxPlot.Visible = 'off';
        vp(ii).WhiskerPlot.Visible = 'off';
        vp(ii).ViolinPlot.LineWidth = 1;
        vp(ii).ScatterPlot.Visible = 'off';
    end
    hold off

    plotter([],[],'k','F','','\rmTime to clearance',[],'days',...
        'AxFontSize',7)
    set(get(gca,'XAxis'),'TickLabelRotation',0)
    set(gca,'XTickLabels',{'WT', [char(916) 'TCR']})
    set(gca,'YScale','linear')
    xlim([0.4 2.6]), ylim([0 20])
    
    % Chronic violin plot
    subplotter(2,4,8,'hratios',hratio2,'lpad',1.1-0.4) %(2,3,6)
    vp = violinplot([ttc_wt_chronic; ttc_dTCR1_chronic],...
        [repmat({'WT1'},[Nsims 1]);...
        repmat({[char(916) 'TCR rep1.']},[Nsims 1])],'GroupOrder',...
        {'WT1',[char(916) 'TCR rep1.']},'ShowData',...
        true);
    vp(1).EdgeColor = wt_col; vp(2).EdgeColor = tdt_col;
    vp(1).ViolinColor = wt_col; vp(2).ViolinColor = tdt_col;
    for ii = 1:2
        vp(ii).MedianPlot.Visible = 'off';
        vp(ii).MeanPlot.Visible = 'on';
        vp(ii).BoxPlot.Visible = 'off';
        vp(ii).WhiskerPlot.Visible = 'off';
        vp(ii).ViolinPlot.LineWidth = 1;
        vp(ii).ScatterPlot.Visible = 'off';
    end
    hold off

    plotter([],[],'k','','','',[],[],'AxFontSize',7)
    set(get(gca,'XAxis'),'TickLabelRotation',0)
    set(gca,'XTickLabels',{'WT', [char(916) 'TCR']})
    set(gca,'YScale','linear','YTick',[0 25 50 75 100 125])
    xlim([0.4 2.6]), ylim([0 125]),
    
    % Compute p-value + plot significance stars
    pvalue_acute = ranksum(ttc_wt_acute,ttc_dTCR1_acute);
    pvalue_chronic = ranksum(ttc_wt_chronic,ttc_dTCR1_chronic); 
    disp("Rank sum p-value (acute):   " + num2str(pvalue_acute))
    disp("Rank sum p-value (chronic): " + num2str(pvalue_chronic))
end

%% Figure S6 - Using a model where kappa_E = from uniform dist

if strcmp(fig_num,'S6')
    w = 17.4;   % width in cm
    h = 15;     % height in cm
    rfk = 5; % reduction factor for ak to reduce size of 3d/surf plots
    rft = 1; % reduction factor for time
    vratios = [0.4 0.3 0.3]; % vertical ratios of subplots
    ms = 4; % marker size
    jitterscale = 1; % in percent of min distance between points

    close(figure(6))
    figurer(6,'w',w,'h',h,'PaperSize','letter')
    
    % Plot specs
    acolor = [137 059 179]./255;
    ccolor = [242 155 041]./255;
    acolor0 = [244 224 255]./255;
    acolor1 = [235 199 255]./255;
    acolor2 = [219 153 255]./255;
    ccolor0 = [255 224 184]./255;
    ccolor1 = [255 225 186]./255;
    ccolor2 = [252 200 129]./255;
    cbar_bounds = [0.7 1];
    cbar_width = 0.01;
    xlimacute = [0 40];
    xlimchronic = [0 120];
    
    % Run simulations for acute and chronic
    deterministic = 'on';
    rP_c = 1.15;
    [ta_wt,sola_wt,pa_wt,ta,sola] = ...
        AvC_2_3('acute','WT',xlimacute,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'deterministic',deterministic,...
        'ICs',[ICa,0.001],'Step',0.5,'Ndet',5e3,'useSeed',true,...
        'kappa_E_min',kappa_E_min,'diststyle_kappa_E','uniform');
    [tc_wt,solc_wt,pc_wt,tc,solc] = ...
        AvC_2_3('chronic','WT',xlimchronic,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',1,'kstyle',kstyle,'deterministic',...
        deterministic,'Ndet',5e3,'useSeed',true,'ICs',[ICc 0.001],...
        'kappa_E_min',kappa_E_min,'diststyle_kappa_E','uniform');

    Va = sola_wt.V; % acute pathogen load at constant step size
    Vc = solc_wt.V; % chronic pathogen load at constant step size
    ka_prop = sola.k_prop;
    kc_prop = solc.k_prop;

    % Compute times to clearance
    ttc_acute = computeTTC(ta_wt,sola_wt.V);
    ttc_chronic = computeTTC(tc_wt,solc_wt.V);
    
    % ------------------------------------------------------------------ %
    % Panel (A) - kappa_E values as function of pMHC reactivity

    % Generate data for the kappa_E, many realizations and average
    % (uncomment to rerun new realizations, data saved for efficiency)
    %{
    N = 25;
    kappas = zeros(N,pa_wt.n);
    for ii = 1:N+1
        if ii == N+1
            [~,~,partemp] =  AvC_2_3('chronic','WT',0.1,...
                'kappa_E',kappa_E,'kappa_E_min',kappa_E_min,...
                'deterministic','on','diststyle_kappa_E','uniform');
            kappa_theoretical = partemp.kappa_E;
            break
        end
        [~,~,partemp] = AvC_2_3('chronic','WT',0.1,...
            'kappa_E',kappa_E,'kappa_E_min',kappa_E_min,...
            'diststyle_kappa_E','uniform');
        kappas(ii,:) = partemp.kappa_E;
    end
    %}
    load('AvC_2_3_kappas_uniform.mat') % saved kappa_E data

    subplotter(3,3,2,'hratios',[0.25 0.50 0.25],'vratios',vratios)
    p1 = plot(pc_wt.ak,kappas','color',[0.7 0.7 0.7],'LineWidth',0.5);
    hold on
    p2 = plotter(pc_wt.ak,kappa_theoretical,'k','A',...
        '\rmpMHC reactivity',['\rmExhaustion rate, ',...
        char(954), '_E'],'','day^{-1}'); hold off
    set(gca,'xscale','log','yscale','linear','xtick',[1e-5 1e-4 1e-3 1e-2])
    xlim([Amin,Amax])
    legend([p1(1),p2],{'Individual','Average'},'fontsize',8,'Location',...
        'southeast')
    legend('boxoff')
    ylim([0 8])

    % ------------------------------------------------------------------ %
    % Panel (B) - acute pathogen load
    subplotter(3,3,4,'vratios',vratios)
    plotter(ta_wt,Va,acolor,'B','\rmTime','Path. load',...
        'days','PFU/mL','AxFontSize',7)
    xlim(xlimacute), hold off
    set(gca,'YScale','log','Ylim',[1 1e6],'YTick',[1 1e2 1e4 1e6],...
        'YMinorTick','on')

    % ------------------------------------------------------------------ %
    % Panel () - acute pMHC reactivity profile
    subplotter(3,3,7,'vratios',vratios)
    plotAffinity(ta,pa_wt,ka_prop,rft,rfk,...
        acolor,acolor0,acolor1,acolor2,cbar_bounds,cbar_width,ttc_acute)
    set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2],'XLim',xlimacute,...
        'YLim',[Amin, Amax])

    % ------------------------------------------------------------------ %
    % Panel (C) - chronic pathogen load
    subplotter(3,3,5,'vratios',vratios)
    plotter(tc_wt,Vc,ccolor,'C','\rmTime','Path. load','days','PFU/mL',...
        'AxFontSize',7)
    xlim(xlimchronic), hold off
    set(gca,'YScale','log','Ylim',[1 1e6],'YTick',[1 1e2 1e4 1e6],...
        'YMinorTick','on')

    % ------------------------------------------------------------------ %
    % Panel () - chronic pMHC reactivity profile
    subplotter(3,3,8,'vratios',vratios)
    plotAffinity(tc,pc_wt,kc_prop,rft,rfk,...
        ccolor,ccolor0,ccolor1,ccolor2,cbar_bounds,cbar_width,ttc_chronic)
    set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2],'XLim',xlimchronic,...
        'YLim',[Amin, Amax])

    % ------------------------------------------------------------------ %
    % Panel (D) - chronic pathogen load
    subplotter(3,3,6,'vratios',vratios)
    % Load the variables (result from AvC_2_3_varyingDists.m)
    effectOfMode = load('AvC_2_3_effectOfMode_uniform.mat');
    Nsims = effectOfMode.N; % number of points over which to vary pars
    n = effectOfMode.n; % number of simulations at each par value
    kmode_var = effectOfMode.kmode_var; % modes
    data_mode = effectOfMode.data_mode; % clearance times per mode

    cmean0 = [158 217 154]; % lightest colour
    cmean1 = [039 115 033]; % darkest colour
    cmeans = zeros(Nsims,3); % matrix to interpolate colours

    % Parse data from cell into matrix
    kmode_ttc = zeros(n,Nsims); % time to clearance, varying the mode
    for ii = 1:Nsims
        cmeans(ii,:) = round((cmean1-cmean0).*(ii-1)./Nsims + cmean0)./255;
        kmode_ttc(:,ii) = data_mode{ii}; % each column = ttc at one mode
        staggerscatter(...
            [nan(1,ii-1) 1./kmode_var(ii) nan(1,Nsims-ii)],...
            [nan(n,ii-1) kmode_ttc(:,ii) nan(n,Nsims-ii)],ms,...
            'markeredgecolor',cmeans(ii,:),...
            'markerfacecolor',cmeans(ii,:),'XScale','log',...
            'YScale','log','XSpread',jitterscale), hold on
    end

    set(gca,'xscale','log','yscale','log')
    plotter([],[],'k','D','\rmpMHC reactivity mode',...
        '\rmTime to clear.',[],'days','AxFontSize',7)
    hold off

    % ------------------------------------------------------------------ %
    % Panel (E) - numerical KO of low-affinity T cells
    % Load the data (gradual removal of T cells with low pMHC-reactivity)
    % (to recompute, run 'AvC_2_3_numericalKO.m' with 'uniform' option 
    % set for 'kappaStyle')
    wt_col = [0.6 0.6 0.6];
    tdt_col = [171, 36, 41]./255;
    numKO = load('AvC_2_3_numericalKO_uniform.mat');
    xcutoff = numKO.xcutoff;
    datattc = numKO.datattc;
    N_KO = numKO.N; % number of cut-off (threshold) values
    n_KO = numKO.n; % number of simulations at each threshold value
    k = numKO.k; % 1/pMHC-reactivity
    
    % Plot
    subplotter(3,3,9,'vratios',vratios)
    colors = zeros(N_KO,3);
    for ii = 1:N_KO
        colors(ii,:) = round((wt_col-tdt_col)*255*(ii-1)./N_KO ...
            + tdt_col.*255)./255;
        staggerscatter([nan(1,ii-1) 1./xcutoff(ii) nan(1,N_KO-ii)],...
            [nan(n_KO,ii-1) datattc{ii} nan(n_KO,N_KO-ii)],4,...
            'markeredgecolor',colors(ii,:),...
            'markerfacecolor',colors(ii,:),'XScale','log',...
            'YScale','log','XSpread',0.5), hold on
        hold on
    end
    set(gca,'xscale','log','yscale','log','ytick',[1e1 1e2 1e3 1e4])
    plotter([],[],'k','E','\rmpMHC reactivity cut-off threshold',...
        '\rmTime to clear.',[],'days','AxFontSize',7)
    fs = get(get(gca,'XLabel'),'FontSize'); % carry forward for 3D plot
    ylim([10 1e4])
    hold off
end

%% Movie S1 -- Time series schema and distributions

if strcmp(fig_num,'M1')

    % Set up plotting properties etc.
    w = 17.8;   % width in cm
    h = 9;     % height in cm
    splty = 2;  % 3 rows
    rfk = 5; % reduction factor for ak to reduce size of 3d/surf plots
    rft = 1; % reduction factor for time

    % Colors
    acolor = [137 059 179]./255;
    ccolor = [242 155 041]./255;
    acolor0 = [244 224 255]./255;
    acolor1 = [235, 199, 255]./255;
    acolor2 = [219, 153, 255]./255;
    ccolor0 = [255 224 184]./255;
    ccolor1 = [255, 225, 186]./255;
    ccolor2 = [252, 200, 129]./255;

    % Viral load time series x-axis limits
    xlimacute = [0 40];
    xlimchronic = [0 120];

    close(figure(6))
    gif_fig = figurer(6,'w',w,'h',h);
    gif_fig.Visible = 'on'; % switch to 'off' for slightly faster time

    % Set parameters and options for running simulations
    deterministic = 'on';

    [ta_wt,sola_wt,pa_wt,ta,sola] = ...
        AvC_2_3('acute','WT',xlimacute(end),'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'deterministic',deterministic,...
        'ICs',[ICa,0.001],'Step',1/3,'kappa_E_min',kappa_E_min);
    [tc_wt,solc_wt,pc_wt,tc,solc] = ...
        AvC_2_3('chronic','WT',xlimchronic(end),'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',1,'kstyle',kstyle,'deterministic',...
        deterministic','ICs',[ICc,0.001],'kappa_E_min',kappa_E_min);

    % Compute times to clearance
    ttc_acute = computeTTC(ta_wt,sola_wt.V);
    ttc_chronic = computeTTC(tc_wt,solc_wt.V);
    
    n = pa_wt.n;
    ak = pa_wt.ak;
    Amin = pa_wt.Amin;
    Amax = pa_wt.Amax;
    k = pa_wt.k;

    % ------------------------------------------------------------------ %
    % (A,B) Plots of acute
    ta_wt = ta_wt; % time vector, variable step size
    t2 = ta; % time vector, constant step size
    Va = sola_wt.V;
    ka_prop = sola.k_prop;

    % ------------------------------------------------------------------ %
    subplotter(splty,3,1,'lpad',1.0+0.3,'rpad',0.65+0.2) % (A), virus
    plotter(ta_wt,Va,acolor), hold on
    ylimacute = get(gca,'YLim');
    plotter([ttc_acute ttc_acute],ylimacute,'k--','A',...
        '\rmTime','','days',[],'AxFontSize',7,'LineWidth',0.5);
    ax1 = gca;
    fs = get(ax1.XLabel,'FontSize');
    ylabel('Pathogen load (PFU/mL)','FontSize',get(ax1.XLabel,'FontSize'))
    set(gca,'XLim',xlimacute,'YScale','log','YLim',[1 1e4],...
        'YMinorTick','off')

    % ------------------------------------------------------------------ %
    subplotter(splty,3,2,'lpad',1.0-0.7,'rpad',0.65+1) % (A), T cell
    [x,y]=meshgrid(t2(rft:rft:end),ak(rfk:rfk:end));
    % z = E(rft:rft:end,rfk:rfk:end)';
    z = zeros(size(x));
    for jj = 1:length(t2(rft:rft:end))
        z(:,jj) = smooth(ka_prop(rft.*jj,rfk:rfk:end)) ...
            ./ smooth(max(ka_prop(rft.*jj,rfk:rfk:end)))*100;
    end

    % cmap_a1 = customcolormap([0 1],[acolor; 1 1 1]);
    cmap_a1 = customcolormap([0 0.5 0.9 1],...
        [acolor2; acolor1; acolor0; [1 1 1]]);
    surf(x,y,z,'edgecolor','none'); view(0,90), hold on
    aka_avg = ka_prop*(1./k);
    plot3(t2,aka_avg,110.*ones(size(t2)),'Color',acolor,...
        'LineWidth',1.5)
    plot3([ttc_acute ttc_acute],[Amin Amax],[110 110],'Color','k',...
        'LineWidth',0.5,'LineStyle','--')
    set(gca,'YScale','log','YTick',[1e-4 1e-2],'Layer','top',...
        'XGrid','off','YGrid','off','XMinorGrid','off',...
        'YMinorGrid','off','Colormap',cmap_a1)
    % xlim(xlimacute)
    plotter([],[],'k','B','\rmTime','\rmpMHC reactivity','days',[],...
        'AxFontSize',7)
    ax2 = gca; ax_pos = get(ax2,'Position');
    cbar = colorbar(ax2,'TickDirection','out','Box','off',...
        'Limits',[0 100],'Ticks',[0 25 50 75 100],...
        'Position',[0.5915 0.6118 0.014 0.3147]);
    ax2.Position = ax_pos;
    ylim([Amin Amax])

    % ------------------------------------------------------------------ %
    % (D,E) Plots of the chronic
    ta_wt = tc_wt;
    t2 = tc;
    Va = solc_wt.V;
    kc_prop = solc.k_prop;

    % ------------------------------------------------------------------ %
    subplotter(splty,3,4,'lpad',1.0+0.3,'rpad',0.65+0.2)
    plotter(ta_wt,Va,ccolor), hold on
    ylimchronic = get(gca,'YLim');
    plotter([ttc_chronic ttc_chronic],ylimchronic,'k--','D',...
        '\rmTime','','days',[],'AxFontSize',7,'LineWidth',0.5);
    ylabel('Pathogen load (PFU/mL)','FontSize',fs)
    ax3 = gca;
    set(ax3,'XLim',xlimchronic,'YScale','log','YLim',[1 1e6],...
        'YMinorTick','off','YTick',[1 1e2 1e4 1e6],'XTick',[0 40 80 120])

    % ------------------------------------------------------------------ %
    subplotter(splty,3,5,'lpad',1.0-0.7,'rpad',0.65+1) % (B), T cell
    [x,y]=meshgrid(t2(rft:rft:end),ak(rfk:rfk:end));
    % z = E(rft:rft:end,rfk:rfk:end)';
    z = zeros(size(x));
    for jj = 1:length(t2(rft:rft:end))
        z(:,jj) = smooth(kc_prop(rft.*jj,rfk:rfk:end)) ...
            ./ smooth(max(kc_prop(rft.*jj,rfk:rfk:end)))*100;
    end

    % cmap_a1 = customcolormap([0 1],[acolor; 1 1 1]);
    cmap_c1 = customcolormap([0 0.5 0.9 1],...
        [ccolor2; ccolor1; ccolor0; [1 1 1]]);
    surf(x,y,z,'edgecolor','none'); view(0,90), hold on
    akc_avg = kc_prop*(1./k);
    plot3(t2,akc_avg,110.*ones(size(t2)),'Color',ccolor,...
        'LineWidth',1.5)
    plot3([ttc_chronic ttc_chronic],[Amin Amax],[110 110],'Color','k',...
        'LineWidth',0.5,'LineStyle','--')
    set(gca,'YScale','log','YTick',[1e-4 1e-2],'Layer','top',...
        'XGrid','off','YGrid','off','XMinorGrid','off',...
        'YMinorGrid','off','Colormap',cmap_c1,'XTick',[0 40 80 120])
    plotter([],[],'k','E','\rmTime','\rmpMHC reactivity','days',[],...
        'AxFontSize',7)
    ax4 = gca; ax_pos = get(ax4,'Position');
    cbar = colorbar(ax4,'TickDirection','out','Box','off',...
        'Limits',[0 100],'Ticks',[0 25 50 75 100],...s
        'Position', [0.5915 0.1118 0.014 0.3147]);
    ax4.Position = ax_pos;
    ylim([Amin Amax])
    xlim([0 t2(end)])

    % Annotations
    annotation('textbox',[0.56 0.928 0.096 0.070],...
        'String',{'% of max'},'LineStyle','none','FontSize',8);
    annotation('textbox',[0.56 0.426 0.096 0.070],...
        'String',{'% of max'},'LineStyle','none','FontSize',8);

    % ------------------------------------------------------------------ %
    % Variable subplots (C,F)
    hold(ax1,'on'), xlim1 = ax1.XLim; ylim1 = ax1.YLim;
    hold(ax2,'on'), xlim2 = ax2.XLim; ylim2 = ax2.YLim;
    hold(ax3,'on'), xlim3 = ax3.XLim; ylim3 = ax3.YLim;
    hold(ax4,'on'), xlim4 = ax4.XLim; ylim4 = ax4.YLim;

    ms = 10;
    timecol = 'r'; % [171, 36, 41]./255;

    for ii = 1:length(ta)
        estimateTimeLeft(ii,1,length(ta),1)

        % Time point in panel A
        l1 = plot(ax1,[ta(ii) ta(ii)],ylim1,'color',timecol);
        p1 = scatter(ax1,ta(ii),sola.V(ii),ms,...
            'MarkerFaceColor',timecol,...
            'MarkerEdgeColor',timecol);

        % Time point in panel B
        l2 = plot3(ax2,[ta(ii) ta(ii)],ylim2,[120 120],'color',timecol);
        p2 = scatter3(ax2,ta(ii),aka_avg(ii),120,ms,...
            'MarkerFaceColor',timecol,...
            'MarkerEdgeColor',timecol);

        % Time point in panel D
        l3 = plot(ax3,[tc(ii) tc(ii)],ylim3,'color',timecol);
        p3 = scatter(ax3,tc(ii),solc.V(ii),ms,...
            'MarkerFaceColor',timecol,...
            'MarkerEdgeColor',timecol);

        % Time point in panel E
        l4 = plot3(ax4,[tc(ii) tc(ii)],ylim4,[120 120],'color',timecol);
        p4 = scatter3(ax4,tc(ii),akc_avg(ii),120,ms,...
            'MarkerFaceColor',timecol,...
            'MarkerEdgeColor',timecol);

        % Acute T cell distribution (C)
        subplotter(splty,3,3,'lpad',1.0-0.1)
        fill([ak' ak(end) ak(1)],...
            [ka_prop(1,:)./max(ka_prop(1,:))*100 0 0],...
            [0.9 0.9 0.9],'LineStyle','none'), hold on
        fill([ak' ak(end) ak(1)],...
            [ka_prop(ii,:)./max(ka_prop(ii,:))*100 0 0],...
            acolor1,'LineStyle','none')
        plotter(ak',ka_prop(ii,:)./max(ka_prop(ii,:))*100,acolor,'C',...
            '\rmpMHC reactivity','',[],[],'AxFontSize',7); hold off
        ylabel('T cell proportion (% max)','FontSize',fs)
        xlim([Amin Amax])
        set(gca,'XScale','log','XTick', [1e-4 1e-2])

        % Chronic T cell distribution (F)
        subplotter(splty,3,6,'lpad',1.0-0.1)
        fill([ak' ak(end) ak(1)],...
            [kc_prop(1,:)./max(kc_prop(1,:))*100 0 0],...
            [0.9 0.9 0.9],'LineStyle','none'), hold on
        fill([ak' ak(end) ak(1)],...
            [kc_prop(ii,:)./max(kc_prop(ii,:))*100 0 0],...
            ccolor1,'LineStyle','none')
        plotter(ak',kc_prop(ii,:)./max(kc_prop(ii,:))*100,ccolor,'F',...
            '\rmpMHC reactivity','','days',[],'AxFontSize',7); hold off
        ylabel('T cell proportion (% max)','FontSize',fs)
        xlim([Amin Amax])
        set(gca,'XScale','log','XTick', [1e-4 1e-2])

        % Annotate the time stamps
        at1 = annotation('textbox',[0.88 0.928 0.15 0.070],...
            'String',{['t = ' num2str(ta(ii), '%.0f') ' days']},...
            'LineStyle','none','FontSize',8,'Color',timecol);
        at2 = annotation('textbox',[0.88 0.426 0.15 0.070],...
            'String',{['t = ' num2str(tc(ii), '%.0f') ' days']},...
            'LineStyle','none','FontSize',8,'Color',timecol);

        % -------------------------------------------------------------- %
        % Export image
        
        exportgraphics(gif_fig,'timeSeries.gif','Append',true,...
            'Resolution',325)

        delete(l1), delete(p1)
        delete(l2), delete(p2)
        delete(l3), delete(p3)
        delete(l4), delete(p4)
        delete(at1), delete(at2)
    end
end

%% Movie S2 -- Nullclines

if strcmp(fig_num,'M2')

    w = 17.8;   % width in cm
    h = 6;    % height in cm

    close(figure(7))
    gif_fig = figurer(7,'w',w,'h',h);
    gif_fig.Visible = 'off';

    % Set parameters and options for running simulations
    deterministic = 'on';
    
    % Define custom time step (need to limit movie size but also need high
    % enough resolution near clearance)
    tc = [0:0.5:20 21:30 32:2:60 60.5:0.5:63.5 ...
        63.55:0.05:64.35 64.36:0.001:64.37 64.38:0.01:64.4 64.5 65:110];

    % Run simulations for chronic
    [~,solc,pc] = ...
        AvC_2_3('chronic','WT',tc,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'deterministic',...
        deterministic,'Step',0.5,'ICs',[ICc,0.001],'kappa_E_min',...
        kappa_E_min,'useSeed',true); % ensure exact same result each time

    % ------------------------------------------------------------------ %
    for ii = 1:length(tc)
        estimateTimeLeft(ii,1,length(tc),1)
        % Plots of the acute
        t = tc(1:ii);
        Vtot = solc.V(1:ii);
        Etot = solc.Etot(1:ii,:);
        ak = pc.ak;

        kappaE_wavg = solc.k_prop(ii,:)*pc.kappa_E;
        k_wavg = solc.k_prop(ii,:)*pc.k;

        % Nullclines
        dVdt = @(V,E) 1./(log(10).*10.^V) .* ...
            ( rP_c .* (10.^V-1) .* (1-(10.^V-1)./Pmax) ...
            - kappa_P .* (10.^E-1) .* ((10.^V-1)./((10.^V-1)+a.*k_wavg)));
        dEdt = @(V,E) 1./(log(10).*10.^E) .* ...
            ( sigma_E ...
            + r_E .* (10.^E-1) .* ((10.^V-1)./((10.^V-1)+k_wavg)) ...
            - delta_E .* (10.^E-1) ...
            - kappaE_wavg .* (10.^E-1) ...
            .* ((10.^V-1)./((10.^V-1)+b.*k_wavg)) ...
            - epsilon .* (10.^E-1) .* (10.^E-1) );

        Vlo = -1;
        Vhi = 7;
        Elo = 1;
        Ehi = 7;

        f = figure('Visible','off');
        Enull = fimplicit(dEdt,[Vlo Vhi Elo Ehi],'MeshDensity',5e2);
        Enull_V = 10.^(Enull.XData)-1;
        Enull_E = 10.^(Enull.YData)-1; hold on

        % Numerically solve for the parabollic V nullcline and parse
        Vnull = fimplicit(dVdt,[Vlo Vhi Elo Ehi],'MeshDensity',5e2);
        Vnull_V = 10.^(Vnull.XData)-1;
        Vnull_E = 10.^(Vnull.YData)-1;
        hold off; close(f)

        load('colours.mat');
        colour = matcol(5,:); % colour choice for trajectory
        ms = 10; % marker size
        figure(7)

        % Log scale
        subplotter(1,2,1)
        plotter(Vnull_V+1,Vnull_E+1,'k'), hold on
        plotter(Enull_V+1,Enull_E+1,[0.6 0.6 0.6])
        scatter(Vtot(end)+1,Etot(end)+1,ms,'MarkerEdgeColor',colour,...
            'MarkerFaceColor',colour)
        p = plotter(Vtot+1,Etot+1,5,'A','\rmPathogen load',...
            '\rm Total eff. T cells','PFU/mL',[]); p.LineWidth = 0.5;
        set(gca,'XScale','log','YScale','log'), hold off
        xlim([10.^Vlo-1 10.^Vhi]), ylim([10.^Elo 10.^Ehi])

        % Linear scale
        subplotter(1,2,2)
        plotter(Vnull_V+1,Vnull_E+1,'k'), hold on
        plotter(Enull_V+1,Enull_E+1,[0.6 0.6 0.6])
        scatter(Vtot(end)+1,Etot(end)+1,ms,'MarkerEdgeColor',colour,...
            'MarkerFaceColor',colour)
        p = plotter(Vtot+1,Etot+1,5,'B','\rmPathogen load',...
            '\rmTotal eff. T cells','PFU/mL',[]); p.LineWidth = 0.5;
        hold off
        xlim([0 1.5e5]), ylim([0 1e6])

        tb = annotation('textbox',...
            'String',"t = " + num2str(tc(ii),'%.1f' + " days"),...
            'Color','r','EdgeColor','none',...
            'Position',[0.36 0.9 0.15 0.05]);
        % -------------------------------------------------------------- %
        % Export image
        exportgraphics(gif_fig,'nullclines.gif','Append',true,...
            'Resolution',325)

        delete(tb)

    end
end

%% Function - compute time to pathogen clearance

function [ttc] = computeTTC(tdata,ydata,tmax)

if nargin < 3
    tmax = tdata(end);
end

tdeplete = tdata(ydata<1);
if isempty(tdeplete)
    ttc = tmax;
    warning(['In computing TTC, found that pathogen not cleared in ' ...
        'alotted time tmax = ' num2str(tmax)])
else
    ttc = tdeplete(1);
end

end

%% Function to plot the pMHC-reactivity

function plotAffinity(t,par,k_prop,rft,rfk,color,color0,color1,color2,...
    cbar_bounds,cbar_width,ttc)

ak = par.ak;
k = par.k;
% k_prop = par.k_prop;

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
if ~isempty(ttc)
    plot3([ttc ttc],[par.Amin par.Amax],[110 110],'Color','k',...
        'LineWidth',0.5,'LineStyle','--')
end
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