function ClusterIndex = getCluster(Wn,Ph,MaxClustNum,MAC_Treshold,PLOT)
% getCluster
%
% INPUTS:
%
% Wn: vector with frequencies.
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
MAC = getMAC(Ph,Ph);  if PLOT==1; plotMAC(MAC,'Antes de Clúster'); end
MAC_linkage = linkage(MAC,'single','euclidean');
ClusterIndex = cluster(MAC_linkage,'MaxClust',MaxClustNum);
ClusterIndex = ClusterIndex==(1:max(ClusterIndex));
%%% -----------------------------------------------------------------------
for ii = 1:size(ClusterIndex,2)
    % Exclude Wn < P25 | Wn > P75
    Location = ClusterIndex(:,ii);
    if sum(Location)>4
        ExcludeThis = false(length(Location),1);
        P25 = prctile(Wn(Location),25);
        P75 = prctile(Wn(Location),75);
        ExcludeThis(Location,1) = Wn(Location)<P25 | Wn(Location)>P75;
        ClusterIndex(logical(ExcludeThis),ii) = false;
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