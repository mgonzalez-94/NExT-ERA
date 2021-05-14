function plotCluster(MOrder,Wn,ClusterIndex,XLabel,YLabel,fYa)
% plotCluster
%
% Response a simple supported beam in 2D.
%
% INPUTS:
%
% MOrder: vector with model order.
% Wn: vector with frequencies.
% ClusterIndex: matrix with logicals, each column is a cluster group.
% XLabel: XLabel string.
% YLabel: XLabel string.
%
% OUTPUTS:
%
% %%%%%%%%%%%%%%%%%%%
% %%% Mateo G. H. %%%
% %%% 2021/05/12  %%%
% %%%%%%%%%%%%%%%%%%%
tic_plotCluster = tic;
%%% -----------------------------------------------------------------------
GraphProp = GraphicalProperties;
%%% -----------------------------------------------------------------------
Fig = figure;
Axes_0 = axes(Fig);
if ~isempty(fYa)
    N = size(fYa,2)-1;
    Clrs = colormap(gray(N+4));
    Clrs([1,2,end-1,end],:) = [];
    Pl1 = plot(Axes_0,fYa(:,1),fYa(:,2:end),...
        'LineStyle','--',...
        'LineWidth',GraphProp.linewidth);
    for ii = 1:N
       Pl1(ii).Color = Clrs(ii,:);
    end
    hold on;
end
Axes_1 = axes(Fig);
N = size(ClusterIndex,2);
Clrs = colormap(hsv(N));
for ii=1:N
    hold on;
    plot(...
        Wn(ClusterIndex(:,ii)),MOrder(ClusterIndex(:,ii)),...
        'LineStyle','none',...
        'Marker','o','MarkerSize',6,...
        'MarkerFaceColor',Clrs(ii,:),'MarkerEdgeColor',0.2*Clrs(ii,:));
    hold on;
end
Axes_1.XLabel.String = XLabel;
Axes_1.YLabel.String = YLabel;
%%% -----------------------------------------------------------------------
if ~isempty(fYa)
    Axes_0.Position = Axes_1.Position;
    Axes_1.Color = 'none';
    %%%
    set(Axes_0,GraphProp.Prop);
    set(Axes_0.XAxis,GraphProp.PropXA);
    set(Axes_0.YAxis,GraphProp.PropYA);
    set(Axes_0.Title,GraphProp.PropT);
    set(Axes_0.XLabel,GraphProp.PropXL);
    set(Axes_0.YLabel,GraphProp.PropYL);
    %%%
    Axes_0.YTickLabels = [];
    Axes_0.XTickLabels = [];
    Axes_0.XLim = Axes_1.XLim;
    Axes_0.XTick = [];
    Axes_0.YTick = [];
    %%%
    set(Fig, 'SizeChangedFcn', {@resize, Axes_0, Axes_1});
end
%%% -----------------------------------------------------------------------
set(Axes_1,GraphProp.Prop);
set(Axes_1.XAxis,GraphProp.PropXA);
set(Axes_1.YAxis,GraphProp.PropYA);
set(Axes_1.Title,GraphProp.PropT);
set(Axes_1.XLabel,GraphProp.PropXL);
set(Axes_1.YLabel,GraphProp.PropYL);
%%% -----------------------------------------------------------------------
disp(['plotCluster: ',num2str(toc(tic_plotCluster),'%.3f')])
end
function resize(~, ~, Axes_0, Axes_1)
Axes_0.Position = Axes_1.Position;
end