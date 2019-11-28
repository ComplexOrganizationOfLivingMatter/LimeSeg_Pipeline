function [cells3dFeatures] = extract3dDescriptors(labelledImage, validCells)
%EXTRACT3DDESCRIPTORS Summary of this function goes here
%   Detailed explanation goes here
cells3dFeatures=[];
for indexCell= 1: max(max(max(labelledImage)))   
oneCell3dFeatures = regionprops3(labelledImage==indexCell, 'PrincipalAxisLength', 'Volume', 'ConvexVolume', 'Solidity', 'SurfaceArea', 'EquivDiameter');
aspectRatio = max(oneCell3dFeatures.PrincipalAxisLength,[],2) ./ min(oneCell3dFeatures.PrincipalAxisLength,[],2);
sphereArea = 4 * pi .* ((oneCell3dFeatures.EquivDiameter) ./ 2) .^ 2;
sphericity = sphereArea ./ oneCell3dFeatures.SurfaceArea;
normalizedVolume = oneCell3dFeatures.Volume/mean(oneCell3dFeatures.Volume);
oneCell3dFeatures = horzcat(oneCell3dFeatures, table(aspectRatio, sphericity, normalizedVolume));
[~,indMax] = max(oneCell3dFeatures.Volume);
cells3dFeatures = [cells3dFeatures; oneCell3dFeatures(indMax,:)];
end
columnIDs = table('Size', size(validCells), 'VariableTypes', {'string'});
columnIDs.Properties.VariableNames = {'ID_Cell'};
columnIDs.ID_Cell = arrayfun(@(x) strcat('cell_', num2str(x)), validCells, 'UniformOutput', false);
cells3dFeatures = horzcat(columnIDs, cells3dFeatures(validCells,:));
end

