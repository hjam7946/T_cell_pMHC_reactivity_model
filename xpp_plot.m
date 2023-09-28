function varargout = xpp_plot(varargin)
% XPP_PLOT    Plot bifurcation diagrams in Matlab that have been saved by
% XPPAUT (Adpated from Mohammad S. Imtiaz, University of Newcastle)
%
%   XPP_PLOT plots the bifurcation diagram with Var1 as the x-axis label
%   and Var2 as the y-axis label (a prompt will require to chose the data
%   file)
%
%   XPP_PLOT(title, xAxisLabel, yAxisLabel) plots the bifurcation diagram
%   with the label specifications (must be character vectors)
%
%   XPP_PLOT(..., xUnit, yUnit) plots the bifurcation diagram with unit
%   specifications in brackets in the axis labels
%
%   opt = XPP_PLOT(...) is a structure containing the line structures
%   plotted in the figure.
%
%   [opt, xData, yData] = XPP_PLOT(...) returns the x-axis and y-axis data
%   values.
%
%   Input parameter specifications: XPP_PLOT(..., 'name', 'value', ...).
%   Let i = {S, U} (denoting stable or unstable) and j = {S, P} (denoting
%   steady-state or periodic orbit, respectively), then available options
%   are:
%
%       'Color_ij'           - the color of the plot lines
%       'Linestyle_ij'       - the plot line styles (e.g. '-', '--', ...)
%       'Linewidth_ij'       - the plot line widths
%       'MarkerSize_ij'      - the size of the markers
%       'MarkerFaceColor_ij' - the inner color of the markers.
%       'FileName'           - the directory for the .dat file, optional
%       'Data'               - the N-by-6 AUTO generated matrix, optional
%
%   To alter the data of the bifurcation diagram by multiplying either the
%   x- or y-data by a factor k, adjust the property 'xFactor' or 'yFactor'
%
%   Example:
%       xpp_plot(..., 'xFactor', 10) multiplies the x-values by 10


p = inputParser;

% Define input parameters (optional ones) -- these rely on the presence
% of bmde_plot() function in the MATLAB path

p.addOptional('title','',@ischar)
p.addOptional('xAxisLabel', 'Var1', @ischar)
p.addOptional('yAxisLabel', 'Var2', @ischar)
p.addOptional('xUnit', [], @(x) ischar(x) || isempty(x) )
p.addOptional('yUnit', [], @(y) ischar(y) || isempty(y) )

p.addParameter('xFactor', 1) % factor to multiply x-axis values
p.addParameter('yFactor', 1) % factor to multiply y-axis values

% Define plot colours
p.addParameter('Color_SS', 'Black') %STABLE STEADY STATE
p.addParameter('Color_US', [0.6 0.6 0.6]) %UNSTABLE STEADY STATE
p.addParameter('Color_SP', [0.4660    0.6740    0.1880]) %STABLE PERIODIC
p.addParameter('Color_UP', 'Blue') %UNSTABLE PERIODIC

% Define plot linestyles
p.addParameter('Linestyle_SS','-');  %STABLE STEADY STATE
p.addParameter('Linestyle_US','--'); %UNSTABLE STEADY STATE
p.addParameter('Linestyle_SP',':');  %STABLE PERIODIC ORBIT
p.addParameter('Linestyle_UP',':');  %UNSTABLE PERIODIC ORBIT

% Define plot linewidths
p.addParameter('Linewidth_SS', 1.5)   %STABLE STEADY STATE
p.addParameter('Linewidth_US', 1)     %UNSTABLE STEADY STATE
p.addParameter('Linewidth_SP', 2)     %STABLE PERIODIC ORBIT
p.addParameter('Linewidth_UP', 1)     %UNSTABLE PERIODIC ORBIT

% Define plot marker sizes
p.addParameter('MarkerSize_SS', 2)    %STABLE STEADY STATE
p.addParameter('MarkerSize_US', 2)    %UNSTABLE STEADY STATE
p.addParameter('MarkerSize_SP', 3)    %STABLE PERIODIC
p.addParameter('MarkerSize_UP', 2)    %UNSTABLE PERIODIC

% Define plot marker face colours
p.addParameter('MarkerFaceColor_SS', 'none')    %STABLE STEADY STATE
p.addParameter('MarkerFaceColor_US', 'none')    %UNSTABLE STEADY STATE
p.addParameter('MarkerFaceColor_SP', 'black')   %STABLE PERIODIC
p.addParameter('MarkerFaceColor_UP', 'none')    %UNSTABLE PERIODIC

% File name, if known
p.addParameter('FileName',[])
p.addParameter('Data',[])

% Parse input parameter values and assign
p.parse(varargin{:});

tit_str = p.Results.title;
Ax_x_title = p.Results.xAxisLabel;
Ax_y_title = p.Results.yAxisLabel;
Ax_x_unit = p.Results.xUnit;
Ax_y_unit = p.Results.yUnit;

xMultFactor = p.Results.xFactor;
yMultFactor = p.Results.yFactor;

C_Ss = p.Results.Color_SS;   
C_Us = p.Results.Color_US;
C_Sp = p.Results.Color_SP;  
C_Up = p.Results.Color_UP;

Lt_Ss = p.Results.Linestyle_SS;
Lt_Us = p.Results.Linestyle_US;
Lt_Sp = p.Results.Linestyle_SP;
Lt_Up = p.Results.Linestyle_UP;

Lw_Ss = p.Results.Linewidth_SS;
Lw_Us = p.Results.Linewidth_US;
Lw_Sp = p.Results.Linewidth_SP;
Lw_Up = p.Results.Linewidth_UP;

Ms_Ss = p.Results.MarkerSize_SS;
Ms_Us = p.Results.MarkerSize_US;
Ms_Sp = p.Results.MarkerSize_SP;
Ms_Up = p.Results.MarkerSize_UP;

Mfc_Ss = p.Results.MarkerFaceColor_SS;
Mfc_Us = p.Results.MarkerFaceColor_US;
Mfc_Sp = p.Results.MarkerFaceColor_SP;
Mfc_Up = p.Results.MarkerFaceColor_UP;

file_name = p.Results.FileName;
data = p.Results.Data;


% Open up the .dat file where the bifurcation data is stored
if isempty(data)
    if isempty(file_name)
        [file_in,path] = uigetfile('*.dat',...
            '.dat file saved by AUTO (XPPAUT) ');
        file_name = [path file_in];
    end
    data = load(file_name); % this is where the data is
else % check that the data fed in is appropriate
    if size(data,2) ~= 6
        error(['Incorrect dimensions in XPP_PLOT, data matrix must be' ...
            'of size N-by-6, i.e. must have size(data,2) == 6'])
    end
end

% if isempty(file_name) && isempty(data)
%     [file_in,path] = uigetfile('*.dat',...
%         '.dat file saved by AUTO (XPPAUT) ');
%     file_name = [path file_in];
% end
% 
% if isempty(data)
%     data = load(file_name); % this is where the data is
% else % check that the data fed in is appropriate
%     if size(data,2) ~= 6
%         error(['Incorrect dimensions in XPP_PLOT, data matrix must be' ...
%             'of size N-by-6, i.e. must have size(data,2) == 6'])
%     end
% end

xData = data(:,1) .* xMultFactor;
yData = zeros(size(data,1),8) * NaN;

for n=1:size(yData,1)
    if n < size(yData,1)
        if data(n,4) == 1 && data(n+1,4) == 1
            yData(n,[1 5]) = data(n,[2 3]);
            
        elseif data(n,4) == 2 && data(n+1,4) == 2
            yData(n,[2 6]) = data(n,[2 3]);
            
        elseif data(n,4) == 3 && data(n+1,4) == 3
            yData(n,[3 7]) = data(n,[2 3]);
            
        elseif data(n,4) == 4 && data(n+1,4) == 4
            yData(n,[4 8]) = data(n,[2 3]);
        end
    else
        if data(n,4) == 1
            yData(n,[1 5]) = data(n,[2 3]);
            
        elseif data(n,4) == 2
            yData(n,[2 6]) = data(n,[2 3]);
            
        elseif data(n,4) == 3
            yData(n,[3 7]) = data(n,[2 3]);
            
        elseif data(n,4) == 4
            yData(n,[4 8]) = data(n,[2 3]);
        end
    end
end

yData = yData .* yMultFactor;

%	-----------------------------------------
h_SS = plot(xData,yData(:,1));
set(h_SS,'color',C_Ss)
set(h_SS,'linestyle',Lt_Ss)
set(h_SS,'linewidth',Lw_Ss)
set(h_SS,'Markersize',Ms_Ss)
set(h_SS,'MarkerFacecolor',Mfc_Ss)
hold on;

%	-----------------------------------------
h_US = plot(xData,yData(:,2));
set(h_US,'color',C_Us)
set(h_US,'linestyle',Lt_Us)
set(h_US,'linewidth',Lw_Us)
set(h_US,'Markersize',Ms_Us)
set(h_US,'MarkerFacecolor',Mfc_Us)

%	-----------------------------------------


h_SP = plot(xData,yData(:,[3 7]));
set(h_SP,'color',C_Sp)
set(h_SP,'linestyle',Lt_Sp)
set(h_SP,'linewidth',Lw_Sp)
set(h_SP,'Markersize',Ms_Sp)
set(h_SP,'MarkerFacecolor',Mfc_Sp)
%	-----------------------------------------


h_UP = plot(xData,yData(:,[4 8]));
set(h_UP,'color',C_Up)
set(h_UP,'linestyle',Lt_Up)
set(h_UP,'linewidth',Lw_Up)
set(h_UP,'Markersize',Ms_Up)
set(h_UP,'MarkerFacecolor',Mfc_Up)
%	-----------------------------------------

grid off;
axis tight

plotter([],[],'k',tit_str,Ax_x_title,Ax_y_title,Ax_x_unit,Ax_y_unit)

% If there is an output argument, return the line structures.
if nargout
    linestructs = struct();
    linestructs.h_SS = h_SS;
    linestructs.h_US = h_US;
    linestructs.h_SP = h_SP;
    linestructs.h_UP = h_UP;
    
    varargout = {linestructs, xData, yData};
end