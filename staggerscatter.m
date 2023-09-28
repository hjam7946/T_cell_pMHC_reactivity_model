function varargout = staggerscatter(x,y,varargin)
% STAGGERSCATTER   Make a violin-plot-esque, categorical-like scatter plot.
% This function uses of the built-in MATLAB function swarmchart.
%
%   STAGGERSCATTER(x,y) will plot a "categorical" scatter plot of the data
%   in the columns of matrix y against the values of x (row or column
%   vector)
%
%   STAGGERSCATTER(x,y,size) will use user-defined marker size for scatter
%   plots
%   
%   STAGGERSCATTER(...,'Name','Value)  allows to use optional input
%   parameters to further specify scatter plot details:
%
%       'XScale':           'linear' (default) or 'log' 
%       'YScale':           'linear' (default) or 'log' 
%       'XSpread':          10 (scatter spread, as percentage of distance
%                           between adjascent x-axis values
%       'MarkerEdgeColor':	'k' (default)
%       'MarkerFaceColor':	'k' (default)

ip = inputParser;

% Define input parameters
ip.addRequired('x')
ip.addRequired('y')
ip.addOptional('Size',5)

ip.addParameter('XScale', 'linear') % or 'log'
ip.addParameter('YScale', 'linear') % or 'log'
ip.addParameter('XSpread', 10) % spread of scatter
ip.addParameter('MarkerEdgeColor','k')
ip.addParameter('MarkerFaceColor','k')


% Parse and assign inputs
ip.parse(x,y,varargin{:});

% N = length(x); % number of "categories" or x-axis values
n = size(y,1); % number of entries per y-axis value

scale = ip.Results.XSpread./100;

x = x(:)'; % maker sure it's a row vector
x = repmat(x(:)',[n,1]);
ms = ip.Results.Size;

if strcmp(ip.Results.XScale,'linear')
    if strcmp(ip.Results.YScale,'linear')
        s = swarmchart(x(:),y(:),ms,...
            'MarkerEdgeColor',ip.Results.MarkerEdgeColor,...
            'MarkerFaceColor',ip.Results.MarkerFaceColor,...
            'XJitterWidth',2.*scale);
    elseif strcmp(ip.Results.YScale,'log')
        ylog = log10(y);
        tempfig = figure('visible','on');
        stemp = swarmchart(x(:),ylog(:),ms,...
            'MarkerEdgeColor',ip.Results.MarkerEdgeColor,...
            'MarkerFaceColor',ip.Results.MarkerFaceColor,...
            'XJitterWidth',2.*scale);
        warning('off','MATLAB:hg:EraseModeIgnored')
        warning('off','MATLAB:structOnObject')
        raw_s = struct(stemp);
        warning('on','MATLAB:hg:EraseModeIgnored')
        warning('on','MATLAB:structOnObject')
        raw_data = raw_s.XYZJittered;
        xnew = raw_data(:,1);
        close(tempfig);
        s = scatter(xnew(:),y(:),ms,...
            'MarkerEdgeColor',ip.Results.MarkerEdgeColor,...
            'MarkerFaceColor',ip.Results.MarkerFaceColor);
        set(gca,'YScale','log')
    end
elseif strcmp(ip.Results.XScale,'log')
    xlog = log10(x);
    tempfig = figure('visible','off');
    if strcmp(ip.Results.YScale,'linear')
        stemp = swarmchart(xlog(:),y(:),ms,...
            'MarkerEdgeColor',ip.Results.MarkerEdgeColor,...
            'MarkerFaceColor',ip.Results.MarkerFaceColor,...
            'XJitterWidth',2.*scale);
    elseif strcmp(ip.Results.YScale,'log')
        ylog = log10(y);
        stemp = swarmchart(xlog(:),ylog(:),ms,...
            'MarkerEdgeColor',ip.Results.MarkerEdgeColor,...
            'MarkerFaceColor',ip.Results.MarkerFaceColor,...
            'XJitterWidth',2.*scale);
    end
    warning('off','MATLAB:hg:EraseModeIgnored')
    warning('off','MATLAB:structOnObject')
    raw_s = struct(stemp);
    warning('on','MATLAB:hg:EraseModeIgnored')
    warning('on','MATLAB:structOnObject')
    raw_data = raw_s.XYZJittered;
    xnew = 10.^(raw_data(:,1));
    close(tempfig);
    s = scatter(xnew(:),y(:),ms,...
        'MarkerEdgeColor',ip.Results.MarkerEdgeColor,...
        'MarkerFaceColor',ip.Results.MarkerFaceColor);
    set(gca,'XScale','log','YScale',ip.Results.YScale)
end

if nargout
    varargout{1} = s;
end

end




