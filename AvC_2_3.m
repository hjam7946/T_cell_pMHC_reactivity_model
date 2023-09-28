function [t, sol, varargout] = AvC_2_3(LCMV, type, varargin)
% AVC_2_3: Generate time-series solutions of the continuum-avidity model of
% LCMV acute versus chronic LCMV (Model version 2.3)
%
%   [t, sol] = AVC_2_3(LCMV,type) -- LCMV indicates which system to solve
%   (enter 'acute' or 'chronic'). For 'type', always use 'WT'.
%
%   Function returns a time vector t, as well as a struct sol containing
%   the solution variables, in order:
%       V       -- t-by-1 vector solution for pathogen load (P in paper)
%       E       -- t-by-n matrix solution for n T cell clones
%       Etot    -- t-by-1 vector solution for total T cell population
%       k_prop  -- t-by-n matrix for proportion of T cell avidities
%       k_avg   -- t-by-1 vector of the average k value in time
%
%   [t, sol] = AVC_2_3(LCMV,type,t_end) -- optional input argument 't_end'
%   determines how long the simulation is run for (default = 100 days). If
%   a vector is used instead of a scalar, then AVC_2_3 returns the solution
%   at the time points indicated by the elements of the vector.
%   
%   [t, sol, param] = AVC_2_3(...) also outputs the parameter structure
%   that stores information about the parameter values (constant or
%   function of avidities) used in the computation
%   
%   [t, sol] = AvC_2_3(...,'Name','Value') allows the user to change
%   individual model parameters or features at a time (for a full list of
%   name/value pairs that can be specified, and their descriptions, see 
%   lines 55-107 of the code below. Some examples:
%
%       [t, sol] = AvC_2_3(...,'delta_E',0.5) - changes the T cell natural
%       turnover rate delta_E
%
%       [t, sol] = AvC_2_3(...,'epsilon',5e-6,'kappa_V',0.2) - changes both
%       rates of T cell inter-cellular competition and pathogen clearance
%
%       [t, sol] = AvC_2_3(...,'deterministic','on','useSeed',true) - 
%       removes variability from the simulations by averaging individual
%       realizations of the exhaustion parameter, kappa_E, using Ndet
%       (default = 5e3) samples, drawn using a rng seed
%
%   [t, sol, param, t_step, sol_step] = AvC_2_3(...) allows the user to
%   also output time/solution vectors at constant step size, indicated by
%   the optional name/value pair 'Step' e.g.
%
%       [t, sol, ~, t_step, sol_step] = AvC_2_3(...,'Step',0.5) will
%       additionally return the solution in constant time intervals of 0.5
%       days (useful when making animations, for example)
%   
%   Created on August 3rd 2020, (C) Hassan (Sam) Jamaleddine 2020

% Define input parser
p = inputParser;

% Define input parameters
p.addRequired('LCMV') % simulating fit to acute or chronic LCMV
p.addRequired('type') % always use 'WT', 'TdT' is a legacy option
p.addOptional('t_end',100) % total integration time (can use vector too)
p.addParameter('Step',1) % step-size for optional t_step/sol_step outputs
p.addParameter('ICs', [0.001 0.001]) % initial conditions for V and Etot
% Note: Pathogen load is scaled by 1e5 (so default is 100 PFU/mL)
% Note: Total effector T cells scaled by 1e6 (so default is 1000 cells),
% this value doesn't matter as much because T cells will be initialized to
% their steady state (S0) in the abscence of pathogen

% Default parameters (note- some names are not accurate/consistent with the
% published version of the paper, since their names were not changed
% throughout the years of model modifications. Refer to names in Table S1
% as well as commented descriptions below for the parameters' actual
% meanings
p.addParameter('n',500); % number of discretized TCR clones
p.addParameter('nmax',500); % number of clones with variable exhaustion
p.addParameter('Amin',4.9860e-05); % minimum avidity
p.addParameter('Amax',3.2911e-02); % maximum avidity
p.addParameter('delta_E',0.2815) % either the value, or the average value
p.addParameter('sigma_E',29.6664) % total value across T cells
p.addParameter('r_E',3.1934) % pathogen induced T cell replication
p.addParameter('kappa_E_min',0.3873) % minimum V-induced exhaustion rate 
p.addParameter('kappa_E',3.3369) % single or average exhaustion rate
p.addParameter('epsilon',3.2395e-06) % inter-cellular T cell competition
p.addParameter('a',0.0018) % scaling factor in pathogen clearance term
p.addParameter('b',38.0927) % scaling factor in exhaustion term
p.addParameter('Vmax',1.1706e+05) % maximum pathogen load (Pmax in paper)
p.addParameter('beta_a',0.7660) % acute pathogen rep. rate (r_P in paper)
p.addParameter('beta_c',1.2201) % chronic pathogen repl. (also r_P)
p.addParameter('kappa_V',0.0494) % pathogen clearance rate by T cells
p.addParameter('omega',1) % weight function for the T cell competition
p.addParameter('psi',1) % relative killing efficacies of Teffs
p.addParameter('k_mean_WT',780.6399) % WT MODE k value (1/reactivity)
p.addParameter('k_mean_TdT',1.7580e+03*2/3) % legacy, do not use
p.addParameter('k_std',1.9141) % span of pMHC-reactivity distribution

p.addParameter('isdist_r_E','off') % constant ('off') or from dist ('on')
p.addParameter('isdist_sigma_E','on') % leave on (multi-clonal model)
p.addParameter('isdist_delta_E','off') % leave off
p.addParameter('isdist_kappa_E','on') % leave on (exhaustion)
p.addParameter('diststyle_kappa_E','exponential') % can do 'uniform' too
p.addParameter('direction','descend') % descend: high for high avidities
p.addParameter('includegif','off') % 'on' if want constant steps in time
p.addParameter('t_end_gif',[]) % end of constant step size simulations
p.addParameter('deterministic','off') % turn on for same result each time
p.addParameter('kstyle','log') % leave as log
p.addParameter('mask',[]) % mask for dTCR repertoire, see Methods
p.addParameter('conserveTot',true) % whether to conserve total T cell #
p.addParameter('Ndet',5e3) % no need to change
p.addParameter('useSeed',false) % use 'true' for identical results
p.addParameter('initialize',true) % whether to initialize T cells w/o P
p.addParameter('initialize_t',5e2) % how long to initialize pre-infection

% Parse the inputs
p.parse(LCMV, type, varargin{:});

% Feed parameters into the param structure to feed into the ODE model
n = p.Results.n; param.n = n; % the number of clones of a given reactivity
nmax = p.Results.nmax; % range of clones with variable exhaustion
param.epsilon = p.Results.epsilon; % effector T cell competition
param.a = p.Results.a; % proportionality constant relating K_V to K_E
param.b = p.Results.b; % proportionality constant relating K_I to K_E
param.Vmax = p.Results.Vmax; % maximum viral load
param.kappa_V = p.Results.kappa_V; % viral clearance rate 0.0260
param.omega = p.Results.omega; % weight function for the T cell competition
param.psi = p.Results.psi; % relative killing efficacies of Teffs

% Check for correct LCMV and assign the value of replication beta
if strcmp(LCMV,'acute')
    param.beta = p.Results.beta_a;
elseif strcmp(LCMV,'chronic')
    param.beta = p.Results.beta_c;
else
    error('Invalid LCMV, enter ''acute'' or ''chronic''')
end

% Check for correct type and assign the location in k for average sigma_E
if strcmp(type,'WT')
    k_mode = p.Results.k_mean_WT;
    k_span = p.Results.k_std;
elseif strcmp(type,'TdT')
    k_mode = p.Results.k_mean_TdT;
    k_span = p.Results.k_std;
else
    error(['Invalid type, enter ''WT'' for wild-type ',...
        'or ''TdT'' for TdT knock-out'])
end

% Define the avidity (and reciprocal avidity) vectors
Amin = p.Results.Amin; % the minimum value of pMHC-reactivity
Amax = p.Results.Amax; % the maximum value of pMHC-reactivity

if strcmp(p.Results.kstyle,'linear')
    ak = logspace(log10(Amax),log10(Amin),n)'; % pMHC-reactivity values
    k = 1./ak; % proportionally reciprocal
elseif strcmp(p.Results.kstyle,'log')
    k = logspace(log10(1./Amax),log10(1./Amin),n)';
    ak = 1./k; % pMHC reactivity
else
    error(['Optional parameter kystle must be ',...
        'either ''linear'' or ''log'''])
end

param.Amin = Amin;
param.Amax = Amax;
param.k = k;
param.ak = ak;

% Assign parameter r_E, the effector T cell expansion rate
if strcmp(p.Results.isdist_r_E,'on')
    dist_r_E = makedist('uniform','lower',1,'upper',p.Results.r_E);
    if strcmp(p.Results.deterministic,'off')
        param.r_E = sort(random(dist_r_E,n,1),p.Results.direction);
    else
        param.r_E = getDeterministic(dist_r_E,1e3,p.Results.direction,n);
    end
elseif strcmp(p.Results.isdist_r_E,'off') % default (should stay that way)
    param.r_E = p.Results.r_E;
else
    error('Optional parameter isdist_r_e must be either ''on'' or ''off''')
end

% Assign parameter sigma_E, the effector T cell source term
if strcmp(p.Results.isdist_sigma_E,'on')
    if strcmp(p.Results.kstyle,'linear')
        dist_sigma_E = makedist('normal','mu',k_mode,'sigma',k_span);
        sigma_E = distribution(dist_sigma_E,p.Results.sigma_E,k,...
            p.Results.mask,p.Results.conserveTot);
        param.sigma_E = sigma_E;
    elseif strcmp(p.Results.kstyle,'log') % default (should stay that way)
        if k_span <= 1 % just a safety net to avoid throwing error
            k_span = 1.5;
            warning('Parameter k_std was changed to be a value > 1!')
        end
        mux = log(k_mode)+(log(k_span)).^2; % k_mode is the desired mode
        sigmax = log(k_span); % k_span is related to the geom. stdev
        dist_sigma_E = makedist('lognormal','mu',mux,'sigma',sigmax);
        sigma_E = distribution(dist_sigma_E,p.Results.sigma_E,k,...
            p.Results.mask,p.Results.conserveTot); % log-normal function
        param.sigma_E = sigma_E;
    end
elseif strcmp(p.Results.isdist_sigma_E,'off')
    param.sigma_E = p.Results.sigma_E ./ n;
else
    error(['Optional parameter isdist_sigma_E must be ',...
        'either ''on'' or ''off'''])
end

% Assign parameter delta_E, the effector T cell turnover rate
if strcmp(p.Results.isdist_delta_E,'on')
    dist_delta_E = makedist('exponential','mu',p.Results.delta_E);
    delta_E_min = p.Results.delta_E_min; % minimum allowable delta_E    
    if strcmp(p.Results.deterministic,'off') 
        param.delta_E = sort(random(dist_delta_E,n,1),...
            p.Results.direction) + delta_E_min; % Teff turnover rate
    else
        param.delta_E = getDeterministic(dist_delta_E,500,...
            p.Results.direction,n) + delta_E_min;
    end
elseif strcmp(p.Results.isdist_delta_E,'off') % default (do not change)
    param.delta_E = p.Results.delta_E;
else
    error(['Optional parameter isdist_delta_E must be ',...
        'either ''on'' or ''off'''])
end


% Assign parameter kappa_E, the effector T cell removal rate by virus
if strcmp(p.Results.isdist_kappa_E,'on')
    % Assign distribution type (exponential, default, or uniform)
    if strcmp(p.Results.diststyle_kappa_E,'exponential') % if exponential
        dist_kappa_E = makedist('exponential','mu',p.Results.kappa_E);
        kappa_E_min = p.Results.kappa_E_min; % minimum allowable kappa_E
    elseif strcmp(p.Results.diststyle_kappa_E,'uniform') % if uniform
        kappa_E_min = p.Results.kappa_E_min; % minimum allowable kappa_E
        dist_kappa_E = makedist('uniform','lower',p.Results.kappa_E_min,...
            'upper',2.*p.Results.kappa_E);
    end

    % Draw nmax parameters from distributions and sort
    % Default: nmax = n = number of unique affinity clones
    if strcmp(p.Results.deterministic,'off')
        param.kappa_E = sort(random(dist_kappa_E,nmax,1),...
            p.Results.direction) + kappa_E_min; % exhaustion rate
    else
        Ndet = p.Results.Ndet;
        param.kappa_E = getDeterministic(dist_kappa_E,Ndet,...
            p.Results.direction,nmax,p.Results.useSeed) + kappa_E_min;
    end

    % If nmax > n (e.g. to expand the range of affinities studied), then
    % let delta = (nmax - n)/2 and set the kappa_E values at the lower
    % bound 1:delta and the upper bound nmax-delta:nmax to be the minimum
    % and maximum values of kappa_E within the range of n clones around the
    % mode, respectively
    if nmax ~= n
        kappa_E_lower = repmat(param.kappa_E(1),[floor((n-nmax)/2) 1]);
        kappa_E_upper = repmat(param.kappa_E(end),[ceil((n-nmax)/2) 1]);
        param.kappa_E = [kappa_E_lower;param.kappa_E;kappa_E_upper];
    end
elseif strcmp(p.Results.isdist_kappa_E,'off')
    param.kappa_E = p.Results.kappa_E;
else
    error(['Optional parameter isdist_kappa_E must be ',...
        'either ''on'' or ''off'''])
end

% Set up initial values of the problem
ICs = p.Results.ICs;

if length(p.Results.ICs)==p.Results.n
    init_V = ICs(1); 
    init_E = ICs(2:end);
else
    if length(p.Results.ICs) < 2
        ICs = repmat(ICs,[1 2]);
    end
    init_V = ICs(1);
    init_Etot = ICs(2);

    % Inital distribution of E
    try
        init_E = distribution(dist_sigma_E,init_Etot,k,...
            p.Results.mask,p.Results.conserveTot)';
    catch
        init_E = init_Etot ./ n .* ones(1,n);
    end
end

init = [init_V init_E];

% Initialize E correctly by solving the ODE without influence of V
if p.Results.initialize
    t_init = p.Results.initialize_t; % how long to initialize for
    try
        [~, sol_init] = ode45(@(t,y)initialize(t,y,param),[0 t_init],...
            init,odeset('NonNegative',1:n+1,'AbsTol',1e-12,'MaxStep',1));
    catch
        keyboard() % for debugging
    end
    init_E = sol_init(end,2:end); % make E at last time point the new IC
    init = [init_V init_E]; % properly rewrite the initial conditions
end

% Solve the ODE system of infection model
t_end = p.Results.t_end;
if length(t_end) == 1
    [t, y] = ode45(@(t,y)AvC_2_3_ode(t,y,param),[0 t_end],init,...
        odeset('NonNegative',1:n+1,'AbsTol',1e-8));%,'MaxStep',1));
elseif length(t_end) > 1
    [t, y] = ode45(@(t,y)AvC_2_3_ode(t,y,param),t_end,init,...
        odeset('NonNegative',1:n+1,'AbsTol',1e-6));%,'MaxStep',1));
else
    error('For ''t_end'', either specify an end time, or a time vector')
end

if nargout > 2
    varargout{1} = param;
end

V = y(:,1) .* 1e5; % pathogen load
E = y(:,2:end) .* 1e6; % number of eff T cells across pMHC-reactivities
Etot = sum(E,2); % sum of effector cells in time
k_prop = E./Etot; % distribution of T cell proportions in time
k_avg = k_prop*k; % matrix multiplication size [t k] * [k 1]

% Spit out the results in the sol structure
sol.V = V;
sol.E = E;
sol.Etot = Etot;
sol.k_prop = k_prop;
sol.k_avg = k_avg;

% If both variable AND constant step size are desired, rerun simulation
% with constant step size
if nargout > 3
    t_step = p.Results.Step; % constant step size
    if ~isempty(p.Results.t_end_gif)
        t_span = 0:t_step:p.Results.t_end_gif;
    else
        t_span = 0:t_step:t_end(end);
    end
    [tgif, solgif] = ode15s(@(t,y)AvC_2_3_ode(t,y,param),t_span,init,...
        odeset('NonNegative',1:n+1,'AbsTol',1e-8,'MaxStep',1));
    
    V = solgif(:,1) .* 1e5; % number of viral particles
    E = solgif(:,2:end) .* 1e6; % number of T cells of all avidities
    Etot = sum(E,2); % sum of effector cells in time
    k_prop = E./Etot; % distribution of T cell proportions in time
    k_avg = k_prop*k; % matrix multiplication size [t k] * [k 1]
    
    % Spit out the results in the sol structure
    solgif = struct();
    solgif.V = V;
    solgif.E = E;
    solgif.Etot = Etot;
    solgif.k_prop = k_prop;
    solgif.k_avg = k_avg;
    
    varargout{2} = tgif;
    varargout{3} = solgif;
end

end

%% System of equations defining the model in resposne to infection

function dy = AvC_2_3_ode(~,y,param)

% Parse parameters
k = param.k;
ak = param.ak;
sigma_E = param.sigma_E;
r_E = param.r_E;
delta_E = param.delta_E;
epsilon = param.epsilon;
a = param.a;
b = param.b;
Vmax = param.Vmax;
beta = param.beta;
kappa_V = param.kappa_V;
kappa_E = param.kappa_E;
omega = param.omega;
psi = param.psi;
n = param.n;

% Make sure solution is treated as a column vector
y = y(:);

% Define the variables
V = y(1)*1e5; % first value is the viral load
E = y(2:end)*1e6; % next n values are the n Teffs of different reactivities

% Define the value k, reciprocal to pMHC reactivity, in n steps
kE = k;
kV = a.*kE;
kI = b.*kE;

% Define the ODEs
dVdt = 1e-5 * ...
    ( beta.*V.*(1-V./Vmax) ...
    - kappa_V .* sum(psi .* E .* (V./(V+kV))) );

dEdt = 1e-6 * ...
    ( sigma_E ...
    + r_E .* E .* (V./(V+kE)) ...
    - delta_E .* E ...
    - kappa_E .* E .* (V./(V+kI)) ...
    - epsilon .* E .* sum(omega .* E) );

dy = [dVdt; dEdt]; % derivative as a function of state variables

end


%% Function that initializes T cells in the ABSENCE of pathogen

function dy = initialize(~,y,param)

% Parse parameters
k = param.k;
ak = param.ak;
sigma_E = param.sigma_E;
r_E = param.r_E;
delta_E = param.delta_E;
epsilon = param.epsilon;
a = param.a;
b = param.b;
Vmax = param.Vmax;
beta = param.beta;
kappa_V = param.kappa_V;
kappa_E = param.kappa_E;
omega = param.omega;
psi = param.psi;
n = param.n;

% Make sure solution is treated as a column vector
y = y(:);

% Define the variables
V = 0; % first value is the pathogen load
E = y(2:end)*1e6; % next n values are the n avidities of Teffs

% Define the avidity reciprocal k in n steps
kE = k;
kV = a.*kE;
kI = b.*kE;

% Define the ODEs (in the ABSENCE of pathogen replication)
dVdt = 0 * ...
    ( beta.*V.*(1-V./Vmax) ...
    - kappa_V .* sum(psi .* E .* (V./(V+kV))) ); % should be 0.

dEdt = 1e-6 * ...
    ( sigma_E ...
    + r_E .* E .* (V./(V+kE)) ...
    - delta_E .* E ...
    - kappa_E .* E .* (V./(V+kI)) ...
    - epsilon .* E .* sum(omega .* E) );

dy = [dVdt; dEdt]; 

end

%% Assign a theoretical distribution linearly

function param = distribution(dist,tot,x,varargin)
% 'dist' is the theoretical distribution from "makedist"
% 'tot' is the probability function distribution normalization factor (i.e.
%   sum of P at all x-axis values = tot)
% 'x' is the x-axis range to define probabilities over
% Can optionally put in variable 'mask' with 0s and 1s to selectively
%   alter parts of the distribution (as done with dTCR rep in Fig. 4C)
%   In this case, 'conserveTot' determines whether to re-scale the
%   probability density in order to conserve the total number of T cells
%   across all pMHC-reactivity values

theoretical = pdf(dist,x); % determine theoretical probability density

% Assign mask below if intend to conserve total number of T cells
if nargin>3
    mask = varargin{1};
    conserveTot = varargin{2};
    if ~isempty(mask) && conserveTot
        theoretical = theoretical.*mask';
    end
end

theoretical = theoretical ./ (sum(theoretical)+eps); % normalize: sum is 1

% Assign mask below if intend to NOT conserve total number of T cells
if nargin>3 && ~isempty(mask) && ~conserveTot
    theoretical = theoretical.*mask';
end

param = tot .* theoretical;
end

%% Obtain deterministic simulation using averaging

% Recall that the simulations have a variable element due to how the
% exhaustion rate is sampled from an exponential distribution then sorted.
% The 'average' behaviour is determined by sampling 'Ndet' times, then 
% taking the average at each pMHC reactivity value.
%
% To produce truly identical simulations (e.g. when doing parameter fitting
% where any stochasticity can lead to difficulties), use a seed by setting
% the name/value pair 'useSeed' to be true

function detPar = getDeterministic(dist,Ndet,spec,n,useSeed)

detPar = zeros(Ndet,n);
if useSeed
    rng(1);
end
for i = 1:Ndet
    detPar(i,:) = sort(random(dist,n,1),spec);
end
detPar = mean(detPar,1)';
if useSeed
    rng shuffle % make sure that a different seed is chosen after the fact
end

end
