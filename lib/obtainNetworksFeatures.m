function [degreeNodesCorrelation,coefCluster,betweennessCentrality] = obtainNetworksFeatures(apicobasal_neighbours,validCells, pathSave)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
addpath(genpath('../MouseMuscleSOD1/Code/Ccs_extraction/Codigo_BCT')); 

%% Networks features
apicobasal_neighbours_binary = zeros(size(apicobasal_neighbours,2),size(apicobasal_neighbours,2));
for nCell=1:size(apicobasal_neighbours,2)
     actualCell=apicobasal_neighbours{1,nCell};
     apicobasal_neighbours_binary(nCell,actualCell)=1; 
end
[~, ~, degreeNodesCorrelation, ~, coefCluster, ~, ~,~,~,~,~,~,~,~,betweennessCentrality]=Prueba_brain(apicobasal_neighbours_binary,apicobasal_neighbours_binary,apicobasal_neighbours_binary,validCells);
save(pathSave, 'coefCluster', 'betweennessCentrality')
end

