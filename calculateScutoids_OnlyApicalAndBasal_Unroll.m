function percentageScutoids=calculateScutoids_OnlyApicalAndBasal_Unroll
%Select path of the salivary gland
selpath=uigetdir('data');
%Load the triangulations of basal and apical surface
load(fullfile(selpath, 'Results', 'unrolledGlands', 'gland_SR_1', 'verticesInfo.mat')); 
newVerticesNeighs2D_basal=newVerticesNeighs2D;
load(fullfile(selpath, 'Results', 'unrolledGlands', 'gland_SR_basal', 'verticesInfo.mat')); 
newVerticesNeighs2D_apical=newVerticesNeighs2D;
load(fullfile(selpath, 'Results', 'valid_cells.mat')); 

%% Calculate neighbours from triangulations
NGs_basal={};
NGs_apical={};
%Apical surface
for k=1:max(max(newVerticesNeighs2D_apical))
  [m,n]=find(newVerticesNeighs2D_apical==k);
  neighs=unique(newVerticesNeighs2D_apical(m,:));
  neighs(neighs==k)=[];
  NGs_apical{k,1}=neighs;
  numberNGs_apical{k,1}=size(NGs_apical{k,1},1);
end
%Basal surface
for k=1:max(max(newVerticesNeighs2D_basal))
  [m,n]=find(newVerticesNeighs2D_basal==k);
  neighs=unique(newVerticesNeighs2D_basal(m,:));
  neighs(neighs==k)=[];
  NGs_basal{k,1}=neighs;
  numberNGs_basal{k,1}=size(NGs_basal{k,1},1);
end
%% Calculate percentage of scutoids in valid cells
 scutoids=cellfun(@(x,y) ~isequal(x,y),NGs_apical,NGs_basal);
 validScutoids=scutoids;
 validScutoids(noValidCells)=[];
 percentageScutoids=sum(validScutoids)/size(validScutoids,1);
end