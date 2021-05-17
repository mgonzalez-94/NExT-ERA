function ClusterIndex = getCluster(varargin)
% getCluster
%
% INPUTS:
%
% Fn: vector with frequencies in Hz.
% Ph: matrix with modal shapes, each row is a coordinate and each column a
%       mode shape.
% MaxClustNum: scalar with max number of clusters groups, usually 4.
% MAC_Treshold: scalar with minimum MAC to remaing in cluster.
% PLOT: plot MAC.
%
% OUTPUTS:
%
% ClusterIndex: matrix with logicals, each column is a cluster group.
%
%
% %%%%%%%%%%%%%%%%%%%
% %%% Mateo G. H. %%%
% %%% 2021/05/12  %%%
% %%%%%%%%%%%%%%%%%%%
tic_getCluster = tic;
%%% -----------------------------------------------------------------------
Fn = varargin{1};
Ph = varargin{2};
MaxClustNum = varargin{3};
MAC_Treshold = varargin{4};
PLOT = varargin{5};
if nargin>5
    FreqAprox = varargin{6};
else
    FreqAprox = [];
end
%%% -----------------------------------------------------------------------
MAC = getMAC(Ph,Ph);  if PLOT==1; plotMAC(MAC,'Antes de Clúster'); end
%%% -----------------------------------------------------------------------
if isempty(FreqAprox)
    MAC_linkage = linkage(MAC,'single','euclidean');
    ClusterIndex = cluster(MAC_linkage,'MaxClust',MaxClustNum);
    ClusterIndex = ClusterIndex==(1:max(ClusterIndex));
else
    ClusterIndex = Fn>=(FreqAprox(:,1)') & Fn<=(FreqAprox(:,2)');
end
%%% -----------------------------------------------------------------------
for ii = 1:size(ClusterIndex,2)
    % Exclude Fn < P25 | Fn > P75
    if isempty(FreqAprox)
        Location = ClusterIndex(:,ii);
        if sum(Location)>4
            ExcludeThis = false(length(Location),1);
            P25 = prctile(Fn(Location),25);
            P75 = prctile(Fn(Location),75);
            ExcludeThis(Location,1) = Fn(Location)<P25 | Fn(Location)>P75;
            ClusterIndex(logical(ExcludeThis),ii) = false;
        end
    end
    % Exclude p50_MAC < MAC_Treshold
    if ~isempty(MAC_Treshold)
        Location = ClusterIndex(:,ii);
        ExcludeThis = false(length(Location),1);
        MAC = getMAC(Ph(:,Location),Ph(:,Location));
        ExcludeThis(Location,1) = prctile(MAC,50,2)<MAC_Treshold;
        ClusterIndex(ExcludeThis,ii) = false;
    end
end
%%% -----------------------------------------------------------------------
disp(['getCluster: ',num2str(toc(tic_getCluster),'%.3f')])
end
