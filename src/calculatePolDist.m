%calculate polygon distribution
function calculatePolDist(H,W,nSeeds,nDiagram,nRealizations)

addpath(genpath('..\Epithelia3D\InSilicoModels\TubularModel\src'))
path2load = ['..\Epithelia3D\InSilicoModels\TubularModel\data\tubularCVT\Data\' num2str(H) 'x' num2str(W) '_' num2str(nSeeds) 'seeds\'];

polygonDis = cell(1,nRealizations);
logNormArea = cell(1,nRealizations);
totalSidesCells = cell(1,nRealizations);
normArea = cell(1,nRealizations);
for nRea = 1:nRealizations
    
    load([path2load 'Image_' num2str(nRea) '_Diagram_' num2str(nDiagram) '.mat'])
    
    [~,sidesCells]=calculateNeighbours(L_original);
    
    noValidCells = unique([L_original(1,:),L_original(end,:)]);
    validCells = setdiff(unique(L_original),noValidCells);
    
    polygonDisImg = calculate_polygon_distribution(sidesCells,validCells);
    polygonDis{nRea} = polygonDisImg(2,:);
    
    
    areaCells = regionprops(L_original,'Area');
    areaCells = cat(1,areaCells.Area);
    areaValidCells = areaCells(validCells);
    logNormArea{nRea} = log10(areaValidCells./(mean(areaValidCells)));
    normArea{nRea} = areaValidCells./(mean(areaValidCells));
    totalSidesCells{nRea} = sidesCells(validCells);

end

polyDist = cell2mat(vertcat(polygonDis{:}));
meanPolyDist = mean(polyDist);
stdPolyDist = std(polyDist);
dispersionLogNormArea = vertcat(logNormArea{:});
dispersionNormArea = vertcat(normArea{:});

relationNormArea_numSides = [horzcat(totalSidesCells{:})',dispersionNormArea];
uniqSides = unique(horzcat(totalSidesCells{:}));
lewis_NormArea = [uniqSides;arrayfun(@(x) mean(relationNormArea_numSides(ismember(relationNormArea_numSides(:,1),x),2)),uniqSides);
    arrayfun(@(x) std(relationNormArea_numSides(ismember(relationNormArea_numSides(:,1),x),2)),uniqSides)];


save([path2load 'polygonDistribution_diag_' num2str(nDiagram) '.mat'],'meanPolyDist','stdPolyDist','polyDist','dispersionLogNormArea','dispersionNormArea','lewis_NormArea')


