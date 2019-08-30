function [cells3dFeatures] = extract3dDescriptors(labelledImage, validCells)
%EXTRACT3DDESCRIPTORS Summary of this function goes here
%   Detailed explanation goes here
cells3dFeatures = regionprops3(labelledImage, 'PrincipalAxisLength', 'Volume', 'ConvexVolume', 'Solidity', 'SurfaceArea', 'EquivDiameter');
aspectRatio = max(cells3dFeatures.PrincipalAxisLength,[],2) ./ min(cells3dFeatures.PrincipalAxisLength,[],2);
sphereArea = 4 * pi .* ((cells3dFeatures.EquivDiameter) ./ 2) .^ 2;
sphericity = sphereArea ./ cells3dFeatures.SurfaceArea;
normalizedVolume = cells3dFeatures.Volume/mean(cells3dFeatures.Volume);

columnIDs = table('Size', size(validCells), 'VariableTypes', {'string'});
columnIDs.Properties.VariableNames = {'ID_Cell'};
columnIDs.ID_Cell = arrayfun(@(x) strcat('cell_', num2str(x)), validCells, 'UniformOutput', false);

cells3dFeatures = horzcat(columnIDs, cells3dFeatures(validCells,:), table(aspectRatio(validCells,:), sphericity(validCells,:), normalizedVolume(validCells), 'VariableNames', {'AspectRatio','Sphericity','normalizedVolume'}));
end

