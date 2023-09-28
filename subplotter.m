% SUBPLOTTER   Make subplots for miscellaneous purposes
%
%   H = SUBPLOTTER(m,n,p), breaks the Figure window into an m-by-n matrix
%   of small axes, selects the p-th axes for the current plot, and returns
%   the axes handle.  The axes are counted along the top row of the Figure
%   window, then the second row, etc.  For example,
%
%         SUBPLOTTER(2,1,1), PLOT(income)
%         SUBPLOTTER(2,1,2), PLOT(outgo)
%
%   plots income on the top half of the window and outgo on the bottom
%   half. If the CurrentAxes is nested in a uipanel the panel is used as
%   the parent for the subplot instead of the current figure.
%
%   This is improved over default version because of how stupid the default
%   version is with respect to labels coming off smaller plots.
%
%   SUBPLOTTER(..., 'Name, 'Value', ...) allows to use optional input
%   parameters to further specify plot details:
%
%       'lpad':     padding, in cm, to the left of each plot (1.1 cm)
%       'rpad':     padding, in cm, to the right of each plot (0.4 cm)
%       'dpad':     padding, in cm, to the bottom of each plot (1.1 cm)
%       'upad':     padding, in cm, to the top of each plot (0.65 cm)
%       'hratios':  fraction of space taken up horizontally by each plot
%       'vratios':  fraction of space taken up vertically by each plot
%
%   Note: for 'hratios' and 'vratios', must specify a vector whose entries
%   sum to 1, and of length corresponding to the number of columns or rows,
%   respectively

function varargout = subplotter(m,n,p,varargin)

ip = inputParser;

% Define input parameters
ip.addRequired('m')
ip.addRequired('n')
ip.addRequired('p')

ip.addParameter('lpad', 1.1, @isnumeric)
ip.addParameter('rpad', 0.4, @isnumeric)
ip.addParameter('dpad', 1.1, @isnumeric)
ip.addParameter('upad', 0.65, @isnumeric)
ip.addParameter('hratios',[])
ip.addParameter('vratios',[])

% Parse and assign inputs
ip.parse(m,n,p,varargin{:});

lpad = ip.Results.lpad; % padding, in cm, to the left of each plot
rpad = ip.Results.rpad; % padding, in cm, to the right of each plot
dpad = ip.Results.dpad; % padding, in cm, to the bottom of each plot
upad = ip.Results.upad; % padding, in cm, to the top of each plot

hratios = ip.Results.hratios; % rel. plot ratios in horizontal direction
vratios = ip.Results.vratios; % rel. plot ratios in vertical direction

% Check that hratios user input is acceptable
if isempty(hratios)
    hratios = ones(1,n)./n; % if nothing specified, assing equal ratios
elseif length(hratios) ~= n || size(hratios(:),2) > 1
    error(['''hratios'' must be a col or row vector of length ' ...
        'equal to the number of columns, ''n'''])
elseif sum(hratios) > 1
    error('Sum of entries in ''hratios'' must be = 1')
end

% Check that vratios user input is acceptable
if isempty(vratios)
    vratios = ones(1,m)./m; % if nothing specified, assign equal ratios
elseif length(vratios) ~= m || size(vratios(:),2) > 1
    error(['''vratios'' must be a col or row vector of length ' ...
        'equal to the number of rows, ''m'''])
elseif sum(vratios) > 1
    error('Sum of entries in ''vratios'' must be = 1')
end

% Check that the n and m are at least 1
if m<1 || n<1
    error('Minimum number of rows and columns is 1')
end

% Check that p is not less than 1 or bigger than nxm
if p<1
    error('Plot index may not be smaller than 1')
elseif p > m*n
    error('Plot index may not exceed the product of rows and columns')
end

fig = gcf;
set(fig,'units','centimeters','paperunits','centimeters')
fig_pos = get(fig,'Position'); % position of the figure with heights in cm

wtot = fig_pos(3); % total width
htot = fig_pos(4); % total height

i = fix((p-1)/n)+1; % which row is it?
j = mod(p-1,n)+1; % which column is it?

plot_width = wtot.*hratios(j); 
plot_height = htot.*vratios(i);

% Obtain normalized (0 to 1) positions and sizes of the subplots
% subplot_xpos = ( (j-1)./n .* total_width + lpad) ./ total_width;
% subplot_ypos = ( (m-i)./m .* total_height + dpad) ./ total_height;
subplot_xpos = ( sum(hratios(1:j-1)) .* wtot + lpad)./ wtot;
subplot_ypos = ( (1-sum(vratios(1:i))) .* htot + dpad)./ htot;
subplot_xlength = ( plot_width - (lpad+rpad)) ./ wtot;
subplot_ylength = ( plot_height - (dpad+upad)) ./ htot;

try
    H = subplot('Position',[subplot_xpos,subplot_ypos,...
        subplot_xlength,subplot_ylength]);
    
    if nargout
        varargout{1} = H;
    end
catch
    warning('Figure is too small to accomodate all these subplots')
end

