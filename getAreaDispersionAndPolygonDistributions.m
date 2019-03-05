clear all
close all

addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files = dir('**/Salivary gland/**/Results/3d_layers_info.mat');

nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
files = files(nonDiscardedFiles);

surfLayers = {'apical','basal'};

polyDistApical = cell(size(files,1),1);
polyDistBasal = cell(size(files,1),1);
logNormAreaApical = cell(size(files,1),1);
logNormAreaBasal = cell(size(files,1),1);
totalSidesCellsApical = cell(size(files,1),1);
totalSidesCellsBasal = cell(size(files,1),1);


for nFile = 1 : size(files,1)
    
    
   
    for nSurf = 1:2
        load([files(nFile).folder '\' surfLayers{nSurf} '\final3DImg.mat'],'img3d')
        load([files(nFile).folder '\' surfLayers{nSurf} '\verticesInfo.mat'],'validCellsFinal','newVerticesNeighs2D')

        neighsCells = cell(1,max(newVerticesNeighs2D(:)));
        for nCell = 1 : max(newVerticesNeighs2D(:))
        	 [nRow,nCol]=find(ismember(newVerticesNeighs2D,nCell));
             cellsNeigh = unique(newVerticesNeighs2D(nRow,:));
             cellsNeigh = cellsNeigh(cellsNeigh~=nCell);
             neighsCells{nCell} = cellsNeigh;             
        end
        sidesCells = cellfun(@(x) length(x), neighsCells);
        [polyDist]=calculate_polygon_distribution(sidesCells,validCellsFinal);
        
        volumePerim = regionprops3(img3d,'Volume');
        areaCells = cat(1,volumePerim.Volume);
        areaValidCells = areaCells(validCellsFinal);
        
        if nSurf == 1
            polyDistApical{nFile} = polyDist(2,:);
            logNormAreaApical{nFile} = log10(areaValidCells./(mean(areaValidCells)));
            totalSidesCellsApical{nFile} = sidesCells(validCellsFinal);
        else
            polyDistBasal{nFile} = polyDist(2,:);
            logNormAreaBasal{nFile} = log10(areaValidCells./(mean(areaValidCells)));
            totalSidesCellsBasal{nFile} = sidesCells(validCellsFinal);
        end
    end
end

polyDistBasal = cell2mat(vertcat(polyDistBasal{:}));
meanPolyDistBasal = mean(polyDistBasal);
stdPolyDistBasal = std(polyDistBasal);

polyDistApical = cell2mat(vertcat(polyDistApical{:}));
meanPolyDistApical = mean(polyDistApical);
stdPolyDistApical = std(polyDistApical);

dispersionNormAreaBasal = vertcat(logNormAreaBasal{:});
dispersionNormAreaApical = vertcat(logNormAreaApical{:});

% figure;histogram(dispersionNormAreaBasal,'Normalization','probability','BinWidth',0.02,'FaceColor',[0,0.4,0.2],'FaceAlpha',1); 
% figure;histogram(dispersionNormAreaApical,'Normalization','probability','BinWidth',0.02,'FaceColor',[0.2,0.8,0],'FaceAlpha',1); 

listPolygons = 3:23;
% figure; bar(listPolygons,meanPolyDistBasal,'FaceColor',[0,0.4,0.2]);  
% hold on; errorbar(listPolygons, meanPolyDistBasal, stdPolyDistBasal, 'k', 'linestyle', 'none'); 
% xlim([3,10])
% ylim([0,0.7])
% 
% figure; bar(listPolygons,meanPolyDistApical,'FaceColor',[0.2,0.8,0]); 
% hold on; errorbar(listPolygons, meanPolyDistApical, stdPolyDistApical, 'k', 'linestyle', 'none'); 
% xlim([3,10])
% ylim([0,0.7])

relationNormArea_numSidesBasal = [horzcat(totalSidesCellsBasal{:})',dispersionNormAreaBasal];
relationNormArea_numSidesApical = [horzcat(totalSidesCellsApical{:})',dispersionNormAreaApical];

uniqSidesBasal = unique(horzcat(totalSidesCellsBasal{:}));
uniqSidesApical = unique(horzcat(totalSidesCellsApical{:}));

lewisBasal_logNormArea = [uniqSidesBasal;arrayfun(@(x) mean(relationNormArea_numSidesBasal(ismember(relationNormArea_numSidesBasal(:,1),x),2)),uniqSidesBasal);
arrayfun(@(x) std(relationNormArea_numSidesBasal(ismember(relationNormArea_numSidesBasal(:,1),x),2)),uniqSidesBasal)];

lewisApical_logNormArea = [uniqSidesApical;arrayfun(@(x) mean(relationNormArea_numSidesApical(ismember(relationNormArea_numSidesApical(:,1),x),2)),uniqSidesApical);
arrayfun(@(x) std(relationNormArea_numSidesApical(ismember(relationNormArea_numSidesApical(:,1),x),2)),uniqSidesApical)];

save(['data\Salivary gland\Wildtype\lewis2D_averagePolygon_AreasDistribution_' date '.mat'],'meanPolyDistBasal','stdPolyDistBasal','meanPolyDistApical','stdPolyDistApical','dispersionNormAreaBasal','dispersionNormAreaApical','listPolygons','lewisBasal_logNormArea','lewisApical_logNormArea')