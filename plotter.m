function varargout = plotter(x,y,varargin)
% PLOTTER   Make plots for miscellaneous purposes
%
%   PLOTTER(x,y) will plot the values of y against x, with no further
%   specifications on the plot (i.e. using MATLAB defaults)
%   
%   PLOTTER(x,y,spec) will plot y against x, using the style specified in
%   input variable spec as follows:
%
%     - If spec is numeric (between 1 and 7), plots with color
%       corresponding to that default color
%     - If spec is numeric (RGB triplet) then RGB color will be plotted
%       (e.g. [0.5 0.5 0.5])
%     - If spec is a character vector, (e.g. 'k-', 'r:'), plots with that
%       specification (must be valid MATLAB syntax for plot style)
%
%   PLOTTER(x,y,spec,lab,xlab,ylab,xunit,yunit) will plot y against x,
%   with the following arguments
%
%     - lab:    plot label, top-left corner by default (e.g. A, B, ...)
%     - xlab:   x-axis label
%     - ylab:   y-axis label
%     - xunit:  unit for x-axis label
%     - yunit:  unit for y-axis label
%
%   plt = PLOTTER(...) returns the line handle for editting the plot(s)
%   after it has been created
%
%   [plt, ax] = PLOTTER(...) returns the line handle for editting as well
%   as the axis handle gca.
%
%   PLOTTER(..., 'Name, 'Value', ...) allows to use optional input
%   parameters to further specify plot details:
%
%       'TickDir':          'out' or 'in' (default)
%       'FontSize':         12 (default)
%       'LineWidth':        1.5 (default)
%       'LineStyle':        '-' (default)
%       'AxLineWidth':      1.5 (default)
%       'AxFontSize':       8 (default)
%       'Interpreter':      'tex' (default)
%       'Box':              'on' (default) or 'off'
%
% See also PLOT, FIGURE

p = inputParser;

% Define input parameters
p.addRequired('x')
p.addRequired('y')
p.addOptional('spec',[],@(s) ischar(s) || isnumeric(s))
p.addOptional('label','', @ischar)
p.addOptional('xAxisLabel', 'Var1', @ischar)
p.addOptional('yAxisLabel', 'Var2', @ischar)
p.addOptional('xUnit', [], @(x) ischar(x) || isempty(x) )
p.addOptional('yUnit', [], @(y) ischar(y) || isempty(y) )

p.addParameter('TickDir', 'out')
p.addParameter('FontSize', 8)
p.addParameter('LineWidth', 1.5)
p.addParameter('LineStyle', '-')
p.addParameter('AxLineWidth', 1)
p.addParameter('AxFontSize', 7)
p.addParameter('AxColor','k')
p.addParameter('TickLength',[0.18 0.18]) % cm, default (norm): [0.01 0.025]
p.addParameter('Interpreter', 'tex')
p.addParameter('Box', 'off')


% Parse and assign inputs
p.parse(x,y,varargin{:});

spec = p.Results.spec;
label = p.Results.label;
xlab = p.Results.xAxisLabel;
ylab = p.Results.yAxisLabel;
xunit = p.Results.xUnit;
yunit = p.Results.yUnit;

tickdir = p.Results.TickDir;
ticklength = p.Results.TickLength; % desired tick length in cm
fs = p.Results.FontSize;
lw = p.Results.LineWidth;
ls = p.Results.LineStyle;
afs = p.Results.AxFontSize;
alw = p.Results.AxLineWidth;
acol = p.Results.AxColor;
interp = p.Results.Interpreter;
boxed = p.Results.Box;

% Set axis handle from gca
ax = gca;

ColorOrder = ...
    [    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

if isempty(spec)
    plt = plot(x,y,'LineWidth',lw,'LineStyle',ls);
elseif ~ischar(spec) && length(spec) == 1
    plt = plot(x,y,'LineWidth',lw,'Color',ColorOrder(spec,:),...
        'LineStyle',ls);
elseif ~ischar(spec) && size(spec,1) == 1 && size(spec,2) == 3
    plt = plot(x,y,'LineWidth',lw,'Color',spec,'LineStyle',ls);
elseif ~ischar(spec) && size(spec,1) > 1 && size(spec,2) == 3
    if size(spec,1) == size(y,2)
        holdoff = ~ishold(ax);
        for ii = 1:size(y,2)
            plot(x,y(:,ii),'LineWidth',lw,'Color',spec(ii,:),...
                'LineStyle',ls), hold on
        end
        if holdoff
            hold off
        end
    else
        error(['If entering a matrix for colour specification into ',...
            'PLOTTER, matrix must have 3 columns and the same ',...
            'number of rows as individual data sets being plotted'])
    end
else
    plt = plot(x,y,spec,'LineWidth',lw,'LineStyle',ls);
end

ax.FontSize = afs;
ax.LineWidth = alw;
ax.TickDir = tickdir;
ax.XColor = acol;
ax.YColor = acol;
ax.ZColor = acol;

% Fix tick length
ax.Units = 'centimeters';
pos = ax.Position;
normticklength = ticklength./max(pos(3:4));
ax.TickLength = normticklength;
ax.Units = 'normalized';

if nargin > 3
    titlestring = ['\bf{',label,'}'];
    
    if isempty(xunit)
        xlabelstring = ['\rm',xlab];
    else
        xlabelstring = ['\rm',xlab, ' (', xunit, ')'];
    end
    
    if isempty(yunit)
        ylabelstring = ['\rm',ylab];
    else
        ylabelstring = ['\rm',ylab, ' (', yunit, ')'];
    end
    
    ax.Units = 'Points';
    ax_height = ax.Position(end); % figure height in pixels
    
    text(-20,ax_height+12,titlestring,'FontSize',fs+2,...
        'Interpreter',interp,'Units','Points')
    xlabel(xlabelstring,'FontSize',fs,'Interpreter',interp)
    ylabel(ylabelstring,'FontSize',fs,'Interpreter',interp)
    
    ax.Units = 'normalized';
end

if strcmp(boxed,'on')
    box on
elseif strcmp(boxed, 'off')
    box off
end

if nargout
    varargout{1} = plt;
    varargout{2} = ax;
end

% set(gcf,'Color',[1 1 1]) % make it white

end