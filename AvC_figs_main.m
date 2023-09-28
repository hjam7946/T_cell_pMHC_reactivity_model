% Run this code to generate the modelling parameters and results shown in
% the figures of the main manuscript (Figures 1B-C, 2A-D, 3, and 4)
%
% Note that, since some results are generated from time-consuming
% simulations, these are instead loaded from .mat files. See where and how
% to regenerate them, on a case-by-case basis, below.
%
% (C) Hassan (Sam) Jamaleddine, 2022

%% Options

fig_num = 5; % choose which figure (with model results) to plot (1-5)

% Load parameters from fit result and simulate
load('AvC - Fits/ga-2022-10-24-at-05-13-49/fitdata.mat');
bestpar=fitspecs.bestpar; % Parameter outputs from genetic algorithm

sigma_E = bestpar(1);       % total thymic input
r_E = bestpar(2);           % T-cell activation/proliferation rate
delta_E = bestpar(3);       % natural t-cell turnover
epsilon = bestpar(4);       % homeostatic deletion rate
Pmax = bestpar(5);          % maximum viral load (Pmax in paper)
rP_c = bestpar(6);          % chronic viral replication rate (rP in paper)
rP_a = bestpar(7);          % acute viral replication rate (also rP)
kappa_E = bestpar(8);       % avg pathogen-dependent T cell exhaustion
kappa_P = bestpar(9);       % pathogen clearance rate by T cells
kmode = 10.^bestpar(10);    % mode value of 1/reactivity (mean of logs)
a = bestpar(11);            % scaling factor in pathogen clearance
b = bestpar(12);            % scaling factor in T cell depletion
kspan = bestpar(13);        % span of the pMHC reactivity distribution
kappa_E_min = bestpar(14);  % minimum exhaustion rate
ICa = bestpar(15);          % initial pathogen load for acute
ICc = bestpar(16);          % initial pathogen load for chronic

Amin = 10.^(-(log10(kmode)+5*log10(kspan))); % define min reactivity
Amax = 10.^(-(log10(kmode)-5*log10(kspan))); % define max reactivity
kstyle = 'log'; % log-normally distributed reactivity (vs. gaussian)
direction = 'descend'; % do not change (or do but it will mess it up)

%% Figure 1 -- Model design (right side, panels B and C)

if fig_num == 1
    w = 8.5;   % width in cm
    h = 10;     % height in cm

    close(figure(1))
    figurer(1,'w',w,'h',h,'PaperSize','letter')
    
    deterministic = 'off';
    
    % Run simulations for acute and chronic, WT
    [ta_wt, sola_wt,pa_wt] = ...
        AvC_2_3('acute','WT',0.1,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'deterministic',deterministic,...
        'ICs',[ICa 0.001],'kappa_E_min',kappa_E_min);
    [tc_wt, solc_wt,pc_wt] = ...
        AvC_2_3('chronic','WT',0.1,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',2,'kstyle',kstyle,'deterministic',...
        deterministic,'ICs',[ICc 0.001],'kappa_E_min',kappa_E_min);

    % Plot sigma_E as a function of ak for both WT and TdTKO
    subplotter(2,1,1)
    fill([pa_wt.ak' pa_wt.ak(end) pa_wt.ak(1)],...
        [pa_wt.sigma_E' 0 0],'r','LineStyle','none'), hold on
    plotter(pa_wt.ak',pa_wt.sigma_E,'k','B',...
        '\rmpMHC reactivity',['\rmThymic input, ',...
        char(963), '_E'],'','cells/day'), hold on
    set(gca,'xscale','log','xtick',[1e-5 1e-4 1e-3 1e-2 1e-1])
    xlim([Amin,Amax])

    % -----------------------------------------------------------------%
    % Generate data for the kappa_E, many realizations and average
    % (uncomment to rerun new realizations, data saved for efficiency)
    %{ 
    N = 25;
    kappas = zeros(N,pa_wt.n);
    for ii = 1:N+1
        if ii == N+1
            [~,~,partemp] =  AvC_2_3('chronic','WT',0.1,...
                'kappa_E',kappa_E,'kappa_E_min',kappa_E_min,...
                'deterministic','on');
            kappa_theoretical = partemp.kappa_E;
            break
        end
        [~,~,partemp] = AvC_2_3('chronic','WT',0.1,...
            'kappa_E',kappa_E,'kappa_E_min',kappa_E_min);
        kappas(ii,:) = partemp.kappa_E;
    end
    %} 
    load('AvC_2_3_kappas.mat') % saved kappas and kappa_theoretical data

     % Plot kappa_E as a function of pMHC-reactivity (colors edited in
     % post-production)
    subplotter(2,1,2)
    p1 = plot(pc_wt.ak,kappas','color',[0.7 0.7 0.7],'LineWidth',0.5);
    hold on
    p2 = plotter(pc_wt.ak,kappa_theoretical,'k','C',...
        '\rmpMHC reactivity',['\rmExhaustion rate, ',...
        char(954), '_E'],'','day^{-1}'); hold off
    set(gca,'xscale','log','yscale','log','xtick',[1e-5 1e-4 1e-3 1e-2])
    xlim([Amin,Amax])
    legend([p1(1),p2],{'Individual','Average'},'fontsize',8,'Location',...
        'southeast')
    legend('boxoff')
    ylim([1e-1 1e2])
end

%% Figure 2 -- Time series schema and distributions

if fig_num == 2
    
    w = 17.4;   % width in cm
    h = 9;     % height in cm
    spltx = 3; % 2 columns
    splty = 2;  % 2 rows
    rfk = 5; % reduction factor for ak to reduce size of 3d/surf plots
    rft = 1; % reduction factor for time
    
    % Colors (acute and chronic, as well as paler versions for heat map)
    acolor = [137 059 179]./255; % acute color (purple)
    ccolor = [242 155 041]./255; % chronic color (orange)
    pcolor = [242 125 041]./255; % persistent color (orange)

    acolor0 = [244 224 255]./255; % palest, acute purple
    acolor1 = [235 199 255]./255; % slightly darker, acute purple
    acolor2 = [219 153 255]./255; % even darker, acute purple
    ccolor0 = [255 224 184]./255; % palest, chronic orange
    ccolor1 = [255 225 186]./255; % slightly darker, chronic orange
    ccolor2 = [252 200 129]./255; % even darker, chronic orange
    pcolor0 = [255 194 184]./255; % palest, persistent orange
    pcolor1 = [255 195 186]./255; % slightly darker, persistent orange
    pcolor2 = [252 170 129]./255; % even darker, persistent orange
    
    % Specs for the colorbar of the heat-map in Figs. 2C and D
    cbar_bounds = [0.7 1]; 
    cbar_width = 0.02;
    
    % How long in time to simulate acute or chronic
    xlimacute = [0 40];
    xlimchronic = [0 100];
    
    close(figure(2))
    figurer(2,'w',w,'h',h,'PaperSize','letter')

    % Run simulations for acute and chronic
    deterministic = 'on'; % ensure using the 'average' exhaustion

    [ta_wt,sola_wt,pa_wt,ta,sola] = ...
        AvC_2_3('acute','WT',xlimacute,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'deterministic',deterministic,...
        'ICs',[ICa,0.001],'Step',0.5,'Ndet',5e3,'useSeed',true,...
        'kappa_E_min',kappa_E_min);
    [tc_wt,solc_wt,pc_wt,tc,solc] = ...
        AvC_2_3('chronic','WT',xlimchronic,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',1,'kstyle',kstyle,'deterministic',...
        deterministic,'Ndet',5e3,'useSeed',true,'ICs',[ICc 0.001],...
        'kappa_E_min',kappa_E_min);
    [tp_wt,solp_wt,pp_wt,tp,solp] = ... % PERSISTENT infection
        AvC_2_3('chronic','WT',xlimchronic,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',1.3,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',1,'kstyle',kstyle,'deterministic',...
        deterministic,'Ndet',5e3,'useSeed',true,'ICs',[ICc 0.001],...
        'kappa_E_min',kappa_E_min);

    % Compute times to clearance
    ttc_acute = computeTTC(ta_wt,sola_wt.V); 
    ttc_chronic = computeTTC(tc_wt,solc_wt.V);

    V = sola_wt.V;
    ka_prop = sola.k_prop;

    n = pa_wt.n;
    ak = pa_wt.ak;
    Amin = pa_wt.Amin;
    Amax = pa_wt.Amax;
    k = pa_wt.k;

    % ------------------------------------------------------------------ %
    subplotter(splty,spltx,1) % (A), acute pathogen load
    plotter(ta_wt,V,acolor), hold on
    ylimacute = get(gca,'YLim');
    plotter([ttc_acute ttc_acute],ylimacute,'k--','A',...
        '\rmTime','','days',[],'AxFontSize',7,'LineWidth',0.5);
    ax1 = gca;
    fs = get(ax1.XLabel,'FontSize');
    ylabel('Path. load (PFU/mL)','FontSize',get(ax1.XLabel,'FontSize'))
    xlim(xlimacute), hold off
    set(ax1,'YScale','log','Ylim',[1 1e6],'YTick',[1 1e2 1e4 1e6],...
        'YMinorTick','off')

    % ------------------------------------------------------------------ %
    subplotter(splty,spltx,4) % (C), acute T cell proportion
    plotAffinity(ta,pa_wt,ka_prop,rft,rfk,acolor,acolor0,acolor1,....
        acolor2,cbar_bounds,cbar_width,[]), hold on
    plot3([ttc_acute ttc_acute],[Amin Amax], [110 110],'k--',...
        'LineWidth',0.5)
    plotter([],[],'k','D','Time','pMHC reactivity','days',[],...
        'AxFontSize',7)
    xlim([0 ta(end)])
    ylim([Amin Amax])
    
    % ------------------------------------------------------------------ %
    % (B) Plots of the chronic 
    V = solc_wt.V;
    kc_prop = solc.k_prop;

    % ------------------------------------------------------------------ %
    subplotter(splty,spltx,2)
    plotter(tc_wt,V,ccolor), hold on
    ylimchronic = get(gca,'YLim');
    plotter([ttc_chronic ttc_chronic],ylimchronic,'k--','B',...
        '\rmTime','','days',[],'AxFontSize',7,'LineWidth',0.5);
    ylabel('Path. load (PFU/mL)','FontSize',fs)
    xlim(xlimchronic), hold off
    set(gca,'XTick',[0 25 50 75 100],'YScale','log','Ylim',[1 1e6],...
        'YTick',[1 1e2 1e4 1e6],'YMinorTick','off')
    
    % ------------------------------------------------------------------ %
    subplotter(splty,spltx,5) % (D), T cell
    plotAffinity(tc,pc_wt,kc_prop,rft,rfk,ccolor,ccolor0,ccolor1,....
        ccolor2,cbar_bounds,cbar_width,[]), hold on
    plot3([ttc_chronic ttc_chronic],[Amin Amax], [110 110],'k--',...
        'LineWidth',0.5)
    plotter([],[],'k','E','Time','pMHC reactivity','days',[],...
        'AxFontSize',7)
    set(gca,'XTick',[0 25 50 75 100])
    xlim([0 tc(end)])
    ylim([Amin Amax])
    
    % ------------------------------------------------------------------ %
    % (C) Plots of the persistent
    V = solp_wt.V;
    kp_prop = solp.k_prop;

    % ------------------------------------------------------------------ %
    subplotter(splty,spltx,3)
    plotter(tp_wt,V,pcolor), hold on
    ylimchronic = get(gca,'YLim');
    plotter([ttc_chronic ttc_chronic],ylimchronic,'k--','C',...
        '\rmTime','','days',[],'AxFontSize',7,'LineWidth',0.5);
    ylabel('Path. load (PFU/mL)','FontSize',fs)
    xlim(xlimchronic), hold off
    set(gca,'XTick',[0 25 50 75 100],'YScale','log','Ylim',[1 1e6],...
        'YTick',[1 1e2 1e4 1e6],'YMinorTick','off')
    
    % ------------------------------------------------------------------ %
    subplotter(splty,spltx,6) % (D), T cell
    plotAffinity(tp,pp_wt,kp_prop,rft,rfk,pcolor,pcolor0,...
        pcolor1,pcolor2,cbar_bounds,cbar_width,[]), hold on
    plotter([],[],'k','F','Time','pMHC reactivity','days',[],...
        'AxFontSize',7)
    set(gca,'XTick',[0 25 50 75 100])
    xlim([0 tc(end)])
    ylim([Amin Amax])
    
    % ------------------------------------------------------------------ %
    
end

%% Figure 3 - Analysis of effect of distribution parameters

if fig_num == 3
    w = 8.5;   % width in cm
    h = 10;    % height in cm
    ms = 4; % marker size for scatter plot
    jitterscale = 1; % in percent of min distance between points
    
    close(figure(3))
    figurer(3,'w',w,'h',h,'PaperSize','letter')

    % Load the variables (result from AvC_2_3_varyingDists.m)
    effectOfDists = load('AvC_2_3_effectOfDists.mat');
    N = effectOfDists.N; % number of points over which to vary pars
    n = effectOfDists.n; % number of simulations at each par value
    kmode_var = effectOfDists.kmode_var; % modes
    kspan_var = effectOfDists.kspan_var; % span
    data_mode = effectOfDists.data_mode; % clearance times per mode
    data_span = effectOfDists.data_span; % clearance times per span
    
    % ------------------------------------------------------------------ %
    % A - Effect of mean
    
    subplotter(2,1,1)
    cmean0 = [158 217 154]; % lightest colour
    cmean1 = [039 115 033]; % darkest colour
    cmeans = zeros(N,3); % matrix to interpolate colours
    
    % Parse data from cell into matrix
    kmode_ttc = zeros(n,N); % time to clearance, varying the mode
    for ii = 1:N
        cmeans(ii,:) = round((cmean1-cmean0).*(ii-1)./N + cmean0)./255;
        kmode_ttc(:,ii) = data_mode{ii}; % each column = ttc at one mode
        staggerscatter(...
            [nan(1,ii-1) 1./kmode_var(ii) nan(1,N-ii)],...
            [nan(n,ii-1) kmode_ttc(:,ii) nan(n,N-ii)],ms,...
            'markeredgecolor',cmeans(ii,:),...
            'markerfacecolor',cmeans(ii,:),'XScale','log',...
            'YScale','log','XSpread',jitterscale), hold on
    end
    
    set(gca,'xscale','log','yscale','log')
    plotter([],[],'k','A','\rmT cell reactivity distribution mode',...
        '\rmTime to clearance',[],'days','AxFontSize',7)
    hold off
    
    % ------------------------------------------------------------------ %
    % B - Effect of span
    subplotter(2,1,2)
    
    cstd0 = [247, 139, 228]; % lightest colour
    cstd1 = [140, 28, 120]; % darkest colour
    cstds = zeros(N,3); % matrix to interpolate colours
    
    % Parse data from cell into matrix
    kspan_ttc = zeros(n,N); % time to clearance, varying the mean/mode
    for ii = 1:N
        cstds(ii,:) = round((cstd1-cstd0).*(ii-1)./N + cstd0)./255;
        kspan_ttc(:,ii) = data_span{ii}; % each column = ttcs at one span
        staggerscatter(...
            [nan(1,ii-1) kspan_var(ii) nan(1,N-ii)],...
            [nan(n,ii-1) kspan_ttc(:,ii) nan(n,N-ii)],ms,...
            'markeredgecolor',cstds(ii,:),...
            'markerfacecolor',cstds(ii,:),...
            'YScale','log','XSpread',jitterscale), hold on
    end
    
    set(gca,'XTick',[1 2 3],'YLim',[10 1e4])
    plotter([],[],'k','B','\rmT cell reactivity distribution spread',...
        '\rmTime to clearance',[],'days','AxFontSize',7)
    hold off
    
    % ------------------------------------------------------------------ %
    % Create inset for the mean/mode (these are schematic only)
    inset1 = axes('Position',[0.700 0.800 0.245 0.141]);
    pdf1 = makedist('normal','mu',-1,'sigma',1);
    pdf2 = makedist('normal','mu',0,'sigma',1);
    pdf3 = makedist('normal','mu',1,'sigma',1);
    fa = 0.8; % face alpha value to use for insets
    
    x = linspace(-5,5,501);
    y1 = pdf(pdf1,x);
    y2 = pdf(pdf2,x);
    y3 = pdf(pdf3,x);
    
    hold on
    fill([x x(end) x(1)], [y1 0 0],cmeans(1,:),'EdgeColor','none',...
        'FaceAlpha',fa)
    plot(x,y1,'k','LineWidth',0.5)
    fill([x x(end) x(1)], [y2 0 0],cmeans(20,:),'EdgeColor','none',...
        'FaceAlpha',fa)
    plot(x,y2,'k','LineWidth',0.5)
    fill([x x(end) x(1)], [y3 0 0],cmeans(end,:),'EdgeColor','none',...
        'FaceAlpha',fa)
    plot(x,y3,'k','LineWidth',0.5)
    hold off
    
    set(inset1,'XTick',[],'YTick',[],'Box','off')
    xlabel(inset1,'\rmpMHC reactivity','FontSize',6)
    ylabel(inset1,'\rmT cell count','FontSize',6)
    
    % ------------------------------------------------------------------ %
    % Create inset for the std
    inset2 = axes('Position',[0.700 0.300 0.245 0.141]);
    pdf1 = makedist('normal','mu',0,'sigma',0.4);
    pdf2 = makedist('normal','mu',0,'sigma',0.8);
    pdf3 = makedist('normal','mu',0,'sigma',1.4);
    
    x = linspace(-5,5,501);
    y1 = pdf(pdf1,x);
    y2 = pdf(pdf2,x);
    y3 = pdf(pdf3,x);
    
    hold on
    fill([x x(end) x(1)], [y1 0 0],cstds(1,:),'EdgeColor','none',...
        'FaceAlpha',fa-0.2)
    plot(x,y1,'k','LineWidth',0.5)
    fill([x x(end) x(1)], [y2 0 0],cstds(15,:),'EdgeColor','none',...
        'FaceAlpha',fa-0.2)
    plot(x,y2,'k','LineWidth',0.5)
    fill([x x(end) x(1)], [y3 0 0],cstds(25,:),'EdgeColor','none',...
        'FaceAlpha',fa-0.2)
    plot(x,y3,'k','LineWidth',0.5)
    hold off
    
    set(inset2,'XTick',[],'YTick',[],'Box','off')
    xlabel(inset2,'\rmpMHC reactivity','FontSize',6)
    ylabel(inset2,'\rmT cell count','FontSize',6)
    
    % Annotations
    annotation('arrow',[0.797 0.860],[0.955 0.955],...
        'HeadSize',4);
    annotation('arrow',[0.839 0.892],[0.412 0.412],...
        'HeadSize',4);
    annotation('arrow',[0.804 0.752],[0.412 0.412],...
        'HeadSize',4);
end

%% Figure 4 - Numerical KO of low-affinity T cells, LCMV-Arm or -Cl13

if fig_num == 4
    w = 11.4;   % width in cm
    h = 9;    % height in cm changed from 9
    tdt_col = [171, 36, 41]./255; % matcol(7,:);
    wt_col = [0.6 0.6 0.6];
    tmax_acute = 30;
    tmax_chronic = 500;
    ms = 3; % marker size, in pts, for scatter plot
    computettcs = false; % change to true if want to compute values fresh
    hratio = [0.45 0.225 0.325];
    
    close(figure(4))
    figurer(4,'w',w,'h',h,'PaperSize','letter')
    % load('colours.mat')
    
    deterministic = 'on';
    
    % Definie deltaTCR rep. mask
    x = logspace(log10(1./Amax),log10(1./Amin),500);
    dTCR_mask = 0.5-0.5*tanh((log10(x)-log10(kmode))/0.25);
    % semilogx(1./x,dTCR_mask) % plot it to check
    
    % Run simulations for acute and chronic, WT
    [ta_wt, sola_wt,pa_wt] = ...
        AvC_2_3('acute','WT',tmax_acute,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'ICs',[ICa,0.001],...
        'deterministic',deterministic,'direction',direction,...
        'kappa_E_min',kappa_E_min);
    [tc_wt, solc_wt,pc_wt] = ...
        AvC_2_3('chronic','WT',tmax_chronic,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',2,'kstyle',kstyle,'deterministic',...
        deterministic,'direction',direction,'ICs',[ICc 0.001],...
        'kappa_E_min',kappa_E_min);
    
    % Run simulations for acute and chronic, dTCR rep.
    [ta_dTCR, sola_dTCR,pa_dTCR] = ...
        AvC_2_3('acute','WT',tmax_acute,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_a',rP_a,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'kstyle',kstyle,'ICs',[ICa,0.001],...
        'deterministic',deterministic,'direction',direction,...
        'mask',dTCR_mask,'kappa_E_min',kappa_E_min);
    [tc_dTCR, solc_dTCR,pc_dTCR] = ...
        AvC_2_3('chronic','WT',tmax_chronic,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',2,'kstyle',kstyle,'deterministic',...
        deterministic,'direction',direction,'mask',dTCR_mask,...
        'ICs',[ICc 0.001],'kappa_E_min',kappa_E_min);
    
    % Run multiple times to get a distribution of times to clearance
    if computettcs == 1
        Nsims = 50; % number of simulations

        % Initialize clearance times
        ttc_wt_acute = zeros(Nsims,1);
        ttc_dTCR_acute = zeros(Nsims,1);
        ttc_wt_chronic = zeros(Nsims,1);
        ttc_dTCR_chronic = zeros(Nsims,1);

        % Add variability specifically to ICs
        ICr = 0.1; % +/- 15% in log space
        ICa = 10.^(log10(ICa).*(1+ICr*randn(Nsims,1)));
        ICc = 10.^(log10(ICc).*(1+ICr*randn(Nsims,1)));

        parfor ii = 1:Nsims
            % estimateTimeLeft(ii,1,Nsims,1)
            
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
            [ta_dTCR_N, sola_dTCR_N, ~] = ...
                AvC_2_3('acute','WT',tmax_acute,'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_a',rP_a,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',2,...
                'kstyle',kstyle,'deterministic','off',...
                'direction',direction,'mask',dTCR_mask,'ICs',...
                [ICa(ii) 0.001],'kappa_E_min',kappa_E_min);
            
            % Chronic
            [tc_wt_N, solc_wt_N,~] = ...
                AvC_2_3('chronic','WT',tmax_chronic,'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_c',rP_c,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',2,...
                'kstyle',kstyle,'deterministic','off',...
                'direction',direction,'ICs',[ICc(ii) 0.001],...
                'kappa_E_min',kappa_E_min);
            [tc_dTCR_N, solc_dTCR_N,~] = ...
                AvC_2_3('chronic','WT',tmax_chronic,'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_c',rP_c,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',2,...
                'kstyle',kstyle,'deterministic','off',...
                'direction',direction,'mask',dTCR_mask,'ICs',...
                [ICc(ii) 0.001],'kappa_E_min',kappa_E_min);
            
            % Times to clearance
            ttc_wt_acute(ii) = ...
                computeTTC(ta_wt_N,sola_wt_N.V,tmax_acute);
            ttc_dTCR_acute(ii) = ...
                computeTTC(ta_dTCR_N,sola_dTCR_N.V,tmax_acute);
            ttc_wt_chronic(ii) = ...
                computeTTC(tc_wt_N,solc_wt_N.V,tmax_chronic);
            ttc_dTCR_chronic(ii) = ...
                computeTTC(tc_dTCR_N,solc_dTCR_N.V,tmax_chronic);
        end
    else
        % Load the data that was computed previously
        % (to recompute, set computettcs = true)
        load('AvC_2_3_WTvsdTCR_lcmv.mat')
    end
    
    % ------------------------------------------------------------------ %
    % B - Numerical KO, data
    
    lpad = 1.1;
    rpad = 0.4;
    dpad = 1.1;
    upad = 0.65;
    
    % Load the data (gradual removal of T cells with low pMHC-reactivity)
    % (to recompute, run 'AvC_2_3_numericalKO.m')
    numKO = load('AvC_2_3_numericalKO.mat');
    xcutoff = numKO.xcutoff;
    datattc = numKO.datattc;
    N_KO = numKO.N; % number of cut-off (threshold) values
    n_KO = numKO.n; % number of simulations at each threshold value
    k = numKO.k; % 1/pMHC-reactivity
    Enaive = numKO.Enaive; % starting T cell population at baseline
    
    %{
    endpt = 20;
    xcutoff = xcutoff([1:end-endpt, end]);
    % par = par([1:end-endpt,end],:);
    datattc = datattc([1:end-endpt,end]);
    N = N-endpt;
    %}
    
    % Plot
    subplotter(2,2,2,'hratio',[hratio(1), hratio(2)+hratio(3)])
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
    plotter([],[],'k','B','\rmpMHC reactivity cut-off threshold',...
        '\rmTime to clear.',[],'days','AxFontSize',7)
    fs = get(get(gca,'XLabel'),'FontSize'); % carry forward for 3D plot
    ylim([10 1e4])
    hold off
    
    % ------------------------------------------------------------------ %
    % A - Numerical KO, baseline ("naive") starting T cell population
    subplotter(2,2,1,'hratio',[hratio(1), hratio(2)+hratio(3)])
    for i = 1:N_KO
        p = plot(1./k, Enaive(i,:));
        p.Color = colors(i,:);
        p.LineWidth = 1.1;
        hold on
        if i == 1
            plotter([],[],'k','A')
        end
    end
    set(gca,'XScale','log','YDir','normal','YTick',[0 1 2]), hold off
    xlim([min(1./k),max(1./k)])
    ylim([0 2])
    set(gca,'XTick',[1e-5 1e-4 1e-3 1e-2],'fontsize',7)
    xlabel('\rmpMHC reactivity','FontSize',fs)
    ylabel('\rmT cell count','FontSize',fs)
    hold off

    % ------------------------------------------------------------------ %
    % C - Distribution, simulation
    subplotter(2,3,4,'hratios',hratio)
    Enaive_wt = solc_wt.E(1,:)./max(solc_wt.E(1,:))*100;
    Enaive_tdt = solc_dTCR.E(1,:)./max(solc_dTCR.E(1,:))*100;
    p2 = plotter(pa_dTCR.ak,Enaive_tdt,tdt_col); hold on
    p1 = plotter(pa_wt.ak,Enaive_wt,wt_col,'C',...
        '\rmpMHC reactivity','\rmT cell count',[],'% of max'); hold off
    set(gca,'xscale','log','XLim',[Amin Amax],...
        'XTick',[1e-5 1e-4 1e-3 1e-2],'FontSize',7)
    l = legend([p1, p2], {'WT', [char(916) 'TCR rep.']},'box','off');
    l.ItemTokenSize = [10 18];
    
    % Inset
    inset = axes('Position',[0.120 0.280 0.080 0.080]);
    plot(pa_dTCR.ak,dTCR_mask,'Color',tdt_col,'LineWidth',1)
    set(inset,'XLim',[Amin Amax],'YLim',[0 1],'XTick',[],'YTick',[0 1],...
        'XScale','log','box','off','FontSize',6)
    xlabel(inset,'pMHC react.','FontSize',6)
    title(inset,{'Select. prob.',...
        '\rm(rel. to WT)'},'FontSize',6)
        
    % ------------------------------------------------------------------ %
    % D,top - Simulation, acute
    subplotter(4,3,8,'hratios',hratio,'dpad',dpad-0.5,'lpad',lpad-0.25,...
        'rpad',rpad-0.25)
    p1 = plotter(ta_dTCR,sola_dTCR.V,tdt_col); hold on
    p2 = plotter(ta_wt,sola_wt.V,wt_col,'D',''); hold off
    set(get(gca,'YLabel'),'String','Pathogen load (PFU/mL)')
    l = legend([p2 p1],{'WT',[char(916) 'TCR rep.']},'box','off');
    l.ItemTokenSize = [10 18];
    ax = gca; ax.FontSize = 7; % ax.YAxis.Exponent = 3;
    set(ax,'XLim',[0 tmax_acute],'YScale','log','YLim',[1 1e4],...
        'YMinorTick','off')
    
    % ------------------------------------------------------------------ %
    % D,bottom - Simulation, chronic
    subplotter(4,3,11,'hratios',hratio,'upad',upad-0.5,'lpad',lpad-0.25,...
        'rpad',rpad-0.25)
    p1 = plotter(tc_dTCR,solc_dTCR.V,tdt_col); hold on
    p2 = plotter(tc_wt,solc_wt.V,wt_col,'','\rmTime',...
        '','days',[]);
    hold off
    set(gca,'XLim',[0 100],'YScale','log','YLim',[1 1e6],...
        'FontSize',7,'YMinorTick','off','YTick',[1 1e2 1e4 1e6])
    
    % ------------------------------------------------------------------ %
    % E - Time to clearance, WT vs. TdT KO (violin plots)
    hratio2 = [hratio 0];
    hratio2(end-1:end) = hratio2(end-1)/2;
    Nsims = length(ttc_wt_chronic);
    addpath Violinplot-Matlab-master/
    
    % Acute violin plot
    subplotter(2,4,7,'hratios',hratio2,'rpad',rpad-0.4) %(2,3,6)
    vp = violinplot([ttc_wt_acute; ttc_dTCR_acute],...
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

    plotter([],[],'k','E','','\rmTime to clearance',[],'days',...
        'AxFontSize',7)
    set(get(gca,'XAxis'),'TickLabelRotation',0)
    set(gca,'XTickLabels',{'WT', [char(916) 'TCR']})
    set(gca,'YScale','linear')
    xlim([0.4 2.6]), ylim([0 20])
    
    % Chronic violin plot
    subplotter(2,4,8,'hratios',hratio2,'lpad',lpad-0.4) %(2,3,6)
    vp = violinplot([ttc_wt_chronic; ttc_dTCR_chronic],...
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
    pvalue_acute = ranksum(ttc_wt_acute,ttc_dTCR_acute); % non-parametric
    pvalue_chronic = ranksum(ttc_wt_chronic,ttc_dTCR_chronic); % non-parametric
    disp("Rank sum p-value (acute):   " + num2str(pvalue_acute))
    disp("Rank sum p-value (chronic): " + num2str(pvalue_chronic))


end

%% Figure 5 - Numerical KO of low-affinity T cells, C. neoformans

if fig_num == 5

    w = 11.4;   % width in cm
    h = 4.5;     % height in cm
    rfk = 5; % reduction factor for ak to reduce size of 3d/surf plots
    rft = 1; % reduction factor for time
    fs = 7; % font size
    hratios = [0.38 0.38 0.24];

    close(figure(9))
    figurer(9,'w',w,'h',h,'PaperSize','letter')
    
    crypto = [4.36e3 1.76e4 1.05e7]; % values digitized from Sionov
    Nsims = 50; % number of simulations
    computeCFUs = false; % change to true if want to compute CFUs fresh

    % Tune parameters to crypto data
    rP_c = 0.4; % replication rate
    Pmax = 1e8; % max pathogen load
    kappa_P = 0.0015; % killing rate by T cells
    ICc = 0.316e3./1e5; % initial pathogen load

    Amin = 10.^(-(log10(kmode)+5*log10(kspan))); % define min reactivity
    Amax = 10.^(-(log10(kmode)-5*log10(kspan))); % define max reactivity
    
    % Define deltaTCR rep. mask as in Fig. 4
    x = logspace(log10(1./Amax),log10(1./Amin),500);
    dTCR_mask = 0.5-0.5*tanh((log10(x)-log10(kmode))/0.15);
    
    % Colors (acute and chronic, as well as paler versions for heat map)
    pcolor = [242 125 041]./255; % persistent color (orange)
    pcolor0 = [255 194 184]./255; % palest, persistent orange
    pcolor1 = [255 195 186]./255; % slightly darker, persistent orange
    pcolor2 = [252 170 129]./255; % even darker, persistent orange
    
    % Specs for the colorbar of the heat-map in Figs. 2C and D
    cbar_bounds = [0.7 1]; 
    cbar_width = 0.04;
    
    % How long in time to simulate acute or chronic
    xlimchronic = [0 30];

    % Run simulations for acute and chronic
    deterministic = 'on'; % ensure using the 'average' exhaustion

    [tc_wt,solc_wt,pc_wt,tc,solc] = ...
        AvC_2_3('chronic','WT',xlimchronic,'sigma_E',sigma_E,...
        'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Pmax,...
        'beta_c',rP_c,'kappa_E',kappa_E,'kappa_V',kappa_P,...
        'k_mean_WT',kmode,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
        'k_std',kspan,'Step',0.25,'kstyle',kstyle,'deterministic',...
        deterministic,'Ndet',5e3,'useSeed',true,'ICs',[ICc 0.001],...
        'kappa_E_min',kappa_E_min,'mask',[]);

    n = pc_wt.n;
    ak = pc_wt.ak;
    Amin = pc_wt.Amin;
    Amax = pc_wt.Amax;
    k = pc_wt.k;
    
    % ------------------------------------------------------------------ %
    % (A) Time series
    V = solc_wt.V;
    
    subplotter(1,3,1,'hratios',hratios)
    scatter([7 14 28], crypto,'o','markeredgecolor',pcolor), hold on
    plotter(tc_wt,V,pcolor,'A','Time','Path. load','days','CFU/lung',...
        'AxFontSize',fs), hold off
    ylimchronic = get(gca,'YLim');
    
    xlim(xlimchronic), %hold off
    set(gca,'XTick',[0 10 20 30],'YScale','log','Ylim',[1e2 1e8],...
        'YTick',[1e2 1e3 1e4 1e5 1e6 1e7 1e8],'YMinorTick','on')
    
    % ------------------------------------------------------------------ %
    % (B) pMHC reactivity profile
    subplotter(1,3,2,'hratios',hratios) 
    kc_prop = solc.k_prop;
    plotAffinity(tc,pc_wt,kc_prop,rft,rfk,pcolor,pcolor0,pcolor1,pcolor2,...
        cbar_bounds,cbar_width,[]), hold on
    plotter([],[],'k','B','Time','pMHC reactivity','days',[],...
        'AxFontSize',fs)
    ylim([Amin Amax])
    xlim([0 tc(end)])
    
    % ------------------------------------------------------------------ %
    % (C) WT vs. dTCR repertoire, CFUs at day
    
    CFU_WT = zeros(Nsims,1);
    CFU_dTCR = zeros(Nsims,1);
    
    if computeCFUs
        for ii = 1:Nsims
            estimateTimeLeft(ii,1,Nsims,1)
            [~,solc_wt] = ...
                AvC_2_3('chronic','WT',[0 28],'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_c',rP_c,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',0.25,...
                'kstyle',kstyle,'deterministic','off','Ndet',5e3,...
                'useSeed',true,'ICs',[ICc 0.001],...
                'kappa_E_min',kappa_E_min);
            [~,solc_dTCR] = ...
                AvC_2_3('chronic','WT',[0 28],'sigma_E',sigma_E,...
                'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,...
                'Vmax',Pmax,'beta_c',rP_c,'kappa_E',kappa_E,...
                'kappa_V',kappa_P,'k_mean_WT',kmode,'a',a,'b',b,...
                'Amin',Amin,'Amax',Amax,'k_std',kspan,'Step',0.25,...
                'kstyle',kstyle,'deterministic','off','Ndet',5e3,...
                'useSeed',true,'ICs',[ICc 0.001],...
                'kappa_E_min',kappa_E_min,'mask',dTCR_mask);
            CFU_WT(ii) = solc_wt.V(end);
            CFU_dTCR(ii) = solc_dTCR.V(end);
        end
    else
        crypto_data = load('AvC_2_3_WTvsdTCR_crypto.mat');
        CFU_WT = crypto_data.CFU_WT;
        CFU_dTCR = crypto_data.CFU_dTCR;
        Nsims = length(CFU_WT);
    end

    subplotter(1,3,3,'hratios',hratios)
    addpath Violinplot-Matlab-master/
    bwidth = 0.4;
    xs = [1 2]; % x-axis values
    tdt_col = [171, 36, 41]./255; % matcol(7,:);
    wt_col = [0.6 0.6 0.6];
    
    % Violin plot
    vp = violinplot([CFU_WT; CFU_dTCR],...
        [repmat({'WT'},[Nsims 1]);...
        repmat({[char(916) 'TCR']},[Nsims 1])],'GroupOrder',...
        {'WT',[char(916) 'TCR']},'ShowData',...
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

    % Compute p-value + plot significance stars
    pvalue = ranksum(CFU_WT,CFU_dTCR); % non-parametric
    % sigstar(xs,pvalue);
    disp("Rank sum p-value (sims): " + num2str(pvalue))

    plotter([],[],'k','C','','\rmPath. load',[],'CFU/lung',...
        'AxFontSize',7)
    % set(gca,'YScale','log')
    set(get(gca,'XAxis'),'TickLabelRotation',0)
    set(gca,'YScale','log','YLim',[6e5 1e7])
    xlim([0.4 2.6])
    %}
end

%% Function - compute time to pathogen clearance

function [ttc] = computeTTC(tdata,ydata,tmax)

if nargin < 3
    tmax = tdata(end);
end

tdeplete = tdata(ydata<1);
if isempty(tdeplete)
    ttc = tmax;
    warning(['In computing time to pathogen clearance, found that ' ...
        'pathogen not cleared in alotted time tmax = ' num2str(tmax)])
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
