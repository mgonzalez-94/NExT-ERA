function plotMAC(varargin)
% plotMAC(MAC)
%
% INPUTS:
%   MAC
%
% OUTPUTS:
%
tic_plotMAC = tic;
%%% -----------------------------------------------------------------------
if nargin>0
    MAC = varargin{1};
end
%%%
Title = '';
if nargin>1
    if ~isempty(varargin{2})
        Title = varargin{2};
    end
end
%%% -----------------------------------------------------------------------
GraphProp = GraphicalProperties;
%%% -----------------------------------------------------------------------
%%% Figura
wid = 13;
hei = 10;
fig = figure('Units','centimeters','Position',[1 1 wid hei]);
fig.Color = [1,1,1];
%%% -----------------------------------------------------------------------
Axes_1 = axes;
Plot_1 = bar3(Axes_1,MAC);
colormap(Axes_1, 'winter');
clb = colorbar;
N   = length(Plot_1);
for kk = 1:N
    Plot_1(kk).CData     = Plot_1(kk).ZData;
    Plot_1(kk).FaceColor = 'interp';
    Plot_1(kk).EdgeColor = 'none';
end
view(-90,90)
Axes_1.YLim = [0.5, size(MAC,1)+0.5];
Axes_1.XLim = [0.5, size(MAC,2)+0.5];
%%% -----------------------------------------------------------------------
set(Axes_1,GraphProp.Prop);
set(Axes_1.XAxis,GraphProp.PropXA);
set(Axes_1.YAxis,GraphProp.PropYA);
set(Axes_1.Title,GraphProp.PropT);
set(Axes_1.XLabel,GraphProp.PropXL);
set(Axes_1.YLabel,GraphProp.PropYL);
%%% -----------------------------------------------------------------------
Axes_1.Box = 'on';
Axes_1.XLabel.String = 'Modo de vibración';
Axes_1.YLabel.String = 'Modo de vibración';
Axes_1.Title.String = Title;
clb.Label.String  = 'MAC';
%%% -----------------------------------------------------------------------
disp(['plotMAC: ',num2str(toc(tic_plotMAC),'%.2f')])
end
%% GraphicalProperties
function GraphProp = GraphicalProperties
GraphProp.fontname            = 'Calibri Light';
GraphProp.fontsize            = 11;
GraphProp.linewidth           = 1.2;
GraphProp.Prop.FontName       = GraphProp.fontname;
GraphProp.Prop.FontSize       = GraphProp.fontsize;
GraphProp.PropYL.FontName     = GraphProp.fontname;
GraphProp.PropYL.FontSize     = GraphProp.fontsize;
GraphProp.PropXL.FontName     = GraphProp.fontname;
GraphProp.PropXL.FontSize     = GraphProp.fontsize;
GraphProp.Prop.GridColor      = [1 1 1]*0.6;
GraphProp.Prop.MinorGridColor = [1 1 1]*0.6;
GraphProp.Prop.XMinorGrid     = 'on';
GraphProp.Prop.YMinorGrid     = 'on';
GraphProp.Prop.XGrid          = 'on';
GraphProp.Prop.YGrid          = 'on';
GraphProp.Prop.Box            = 'on';
GraphProp.PropYA.Color        = [1 1 1]*0;
GraphProp.PropXA.Color        = [1 1 1]*0;
GraphProp.PropT.FontName      = GraphProp.fontname;
GraphProp.PropT.FontSize      = GraphProp.fontsize;
% set(Axes_1,GraphProp.Prop);
% set(Axes_1.XAxis,GraphProp.PropXA);
% set(Axes_1.YAxis,GraphProp.PropYA);
% set(Axes_1.Title,GraphProp.PropT);
% set(Axes_1.XLabel,GraphProp.PropXL);
% set(Axes_1.YLabel,GraphProp.PropYL);
end