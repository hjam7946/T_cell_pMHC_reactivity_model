% Code used exclusively by AvC_2_3_fitting which uses genetic algorithm to
% fit parameters. This code generates the model output at the set of
% parameters being sampled by the genetic algorithm.

function [t,V,par] = AvC_2_3_getModelPrediction(time_data,sigma_E,r_E, ... 
    delta_E,epsilon, Vmax, beta, kappa_E, kappa_V, K_E, a, b, stdev,...
    kappa_E_min, ICa, ICc, type)

ICs = [0.0001 0.001];
if strcmp(type,'acute')
    ICs = [ICa 0.001];
elseif strcmp(type,'chronic')
    ICs = [ICc 0.001];
end

Amin = 10.^(-(log10(K_E)+5*log10(stdev)));
Amax = 10.^(-(log10(K_E)-5*log10(stdev)));
[t, sol, par] = AvC_2_3('chronic','WT',time_data,'sigma_E',sigma_E,...
    'r_E',r_E,'delta_E',delta_E,'epsilon',epsilon,'Vmax',Vmax,...
    'beta_c',beta,'kappa_E',kappa_E,'kappa_V',kappa_V,...
    'k_mean_WT',K_E,'a',a,'b',b,'Amin',Amin,'Amax',Amax,...
    'deterministic','on','ICs',ICs,'k_std',stdev,'kstyle','log',...
    'kappa_E_min', kappa_E_min,'useSeed',true);
V = sol.V;
par.Ess = sol.E(1,:);
par.k_prop = sol.k_prop;
end
