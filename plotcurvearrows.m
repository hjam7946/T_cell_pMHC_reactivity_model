function varargout = plotcurvearrows(xdata,ydata,varargin)
% PLOTCURVEARROWS   Plot arrows on top of curves to indicate direction
%
%   PLOTCURVEARROWS(x,y) will plot 5 roughly evenly spaced arrows along the
%   curve with x and y data, exluding endpoints unless otherwise specified.
%   This uses an arclength calculation to determine spacing.
%   
%   PLOTCURVEARROWS(x,y,N) will plot N arrows along the curve
%
%   arrows = PLOTCURVEARROWS(...) returns a cell of size N-by-1 with
%   entries corresponding to the arrow objects
%
%   PLOTCURVEARROWS(...,'Name','Value') allows for further optional
%   parameter inputs.
%
%       'InclEndpoints':    'yes' or 'no' (default)
%       'HeadSize':         5 (default)
%       'HeadStyle':        'cback3' (refer to other MATLAB styles)
%       'Color':            'black'

ip = inputParser;

ip.addRequired('xdata')
ip.addRequired('ydata')
ip.addOptional('N',5,@isnumeric)
ip.addParameter('InclEndpoints','no',@ischar)
ip.addParameter('HeadSize',5,@isnumeric)
ip.addParameter('HeadStyle','cback3')
ip.addParameter('Color','k')

ip.parse(xdata,ydata,varargin{:})

numArrows = ip.Results.N;
if strcmp(ip.Results.InclEndpoints,'yes')
    inclEndPts = true;
elseif strcmp(ip.Results.InclEndpoints,'no')
    inclEndPts = false;
else
    error('Input ''yes'' or ''no'' for ''InclEndpoints'' parameter')
end
size = ip.Results.HeadSize;
style = ip.Results.HeadStyle;
color = ip.Results.Color;

% Check that data is long enough and of same size
if length(xdata) ~= length(ydata)
    error('Data must be of the same length')
elseif length(xdata) < numArrows+1
    error('More errors than the data allows for')
end

% Compute the arclength
arcLengths = computeArcLength(xdata,ydata);
totArcLength = arcLengths(end);

% Index where to put the arrows based on the number
if ~inclEndPts
    targetALs = linspace(totArcLength./(numArrows+1),...
        totArcLength.*(1-1./(numArrows+1)),numArrows);
elseif inclEndPts
    targetALs = linspace(0,totArcLength,numArrows);
end

arrows = cell(numArrows,1);

for ii = 1:numArrows
    index = find(abs(arcLengths-targetALs(ii)) == ...
        min(abs(arcLengths-targetALs(ii))));
    arrow = annotation('arrow');
    arrow.Position = getPosition(xdata(index), ydata(index), ...
        diff(xdata(index:index+1)), diff(ydata(index:index+1)));
    arrow.HeadSize = size;
    arrow.HeadStyle = style;
    arrow.Color = color;
    
    arrows{ii} = arrow;
end

if nargout
    varargout{1} = arrows;
end

end

%% Function to compute the arclength vector, N-1 in length

function arcLength = computeArcLength(xdata,ydata)
dx = diff(xdata);
dy = diff(ydata);
arcLength = cumsum(sqrt(dx.^2 + dy.^2));
end

%% Get the figure panel position values corresponding to axis position

function position = getPosition(x,y,diffx,diffy)

ax = gca;

origin = [x,y];
origin = origin - [ax.XLim(1),ax.YLim(1)];
origin = origin ./ [diff(ax.XLim),diff(ax.YLim)];
origin = origin .* ax.Position(3:4) + ax.Position(1:2);

direction = [diffx,diffy];
direction = direction ./ [diff(ax.XLim),diff(ax.YLim)];
direction = direction .* ax.Position(3:4);

position = [origin, direction];

end