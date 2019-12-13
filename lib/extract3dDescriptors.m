function [cells3dFeatures] = extract3dDescriptors(labelledImage, validCells)
%EXTRACT3DDESCRIPTORS Summary of this function goes here
%   Detailed explanation goes here
cells3dFeatures=[];
for indexCell= 1: max(max(max(labelledImage)))   
    oneCell3dFeatures = regionprops3(labelledImage==indexCell, 'PrincipalAxisLength', 'Volume', 'ConvexVolume', 'Solidity', 'SurfaceArea', 'EquivDiameter');
    
    indMax = 1;
    if size(oneCell3dFeatures, 1) > 1
        [~,indMax] = max(oneCell3dFeatures.Volume);
        oneCell3dFeatures = oneCell3dFeatures(indMax,:);
    end
    actualImg = bwlabeln(labelledImage==indexCell);
    [x, y, z] = ind2sub(size(labelledImage), find(actualImg==indMax));
    DT = delaunayTriangulation(x, y, z);
    [~, convexVolume] = convexHull(DT);
    oneCell3dFeatures.ConvexVolume = convexVolume;
    oneCell3dFeatures.Solidity = sum(actualImg(:)==indMax) / convexVolume;
    aspectRatio = max(oneCell3dFeatures.PrincipalAxisLength,[],2) ./ min(oneCell3dFeatures.PrincipalAxisLength,[],2);
    sphereSurfaceArea = pi^(1/3) * (6*oneCell3dFeatures.Volume)^(2/3); 
    sphericity = sphereSurfaceArea ./ oneCell3dFeatures.SurfaceArea;
    normalizedVolume = oneCell3dFeatures.Volume;
    oneCell3dFeatures.EquivDiameter = sphereSurfaceArea;
    cells3dFeatures = [cells3dFeatures; horzcat(oneCell3dFeatures, table(aspectRatio, sphericity, normalizedVolume))];
end
cells3dFeatures.normalizedVolume = arrayfun(@(x) x/mean(cells3dFeatures.Volume), cells3dFeatures.normalizedVolume);

columnIDs = table('Size', size(validCells), 'VariableTypes', {'string'});
columnIDs.Properties.VariableNames = {'ID_Cell'};
columnIDs.ID_Cell = arrayfun(@(x) strcat('cell_', num2str(x)), validCells, 'UniformOutput', false);
cells3dFeatures = horzcat(columnIDs, cells3dFeatures(validCells,:));
end

