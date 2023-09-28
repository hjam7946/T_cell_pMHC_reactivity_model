function varargout = figurer(num,varargin)
% FIGURER(num)     Generate a MATLAB figure, using the built-in function
% FIGURE but with easier width and height specifications. Syntax:
%
%    FIGURER(num,'Name','Value') generates figure number num, with optional
%    parameters being the 'width' (default = 17.8 cm), 'height' (default =
%    6 cm) and 'Position' (default = [0cm 20cm] for top-left corner)


ip = inputParser;

ip.addRequired('num',@isnumeric)
ip.addParameter('Width',17.8,@isnumeric)
ip.addParameter('Height',6,@isnumeric)
ip.addParameter('Position',[],@isnumeric)
ip.addParameter('PaperSize','plot')
ip.addParameter('Color','w')

ip.parse(num, varargin{:});
w = ip.Results.Width;
h = ip.Results.Height;
p = ip.Results.Position;
color = ip.Results.Color;
papersize = ip.Results.PaperSize;

if isempty(p)
    set(0,'units','centimeters')
    screensize = get(0,'screensize');
    p = [0 screensize(4)-h-2.25];
end

f = figure(num);
if strcmp(papersize,'plot')
    papersize = [w+eps h+eps];
elseif strcmp(papersize,'letter')
    papersize = [8.5*2.54 11*2.54];
else
    if size(papersize(:),2) > 1 || size(papersize(:),1) ~= 2
        error(['Variable ''PaperSize'' must be set to ''plot'' '...
            '(to be the same size as the figure size), ''letter'' '...
            '(to be the size of a letter page, 8.5''''-by-11'''', '...
            'or to a 1-by-2 vector specifying custom length in cm.'])
    end
end
    
set(f,'units','centimeters','position',[p(1) p(2) w h],...
    'paperunits','centimeters','papersize',papersize,'color',color)

if nargout
    varargout{1} = f;
end