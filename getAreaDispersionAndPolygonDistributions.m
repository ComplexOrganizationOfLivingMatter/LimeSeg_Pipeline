clear all
close all

% addpath(genpath('src'))
% addpath(genpath('lib'))
% addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));
% 
% files = dir('**/Salivary gland/**/Results/3d_layers_info.mat');
% 
% nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
% files = files(nonDiscardedFiles);
% 
% surfLayers = {'apical','basal'};
% 
% polyDistApical = cell(size(files,1),1);
% polyDistBasal = cell(size(files,1),1);
% logNormAreaApical = cell(size(files,1),1);
% logNormAreaBasal = cell(size(files,1),1);
% normAreaApical = cell(size(files,1),1);
% normAreaBasal = cell(size(files,1),1);
% totalSidesCellsApical = cell(size(files,1),1);
% totalSidesCellsBasal = cell(size(files,1),1);
% 
% 
% for nFile = 1 : size(files,1)
%     
%     
%    
%     for nSurf = 1:2
%         load([files(nFile).folder '\' surfLayers{nSurf} '\final3DImg.mat'],'img3d')
%         load([files(nFile).folder '\' surfLayers{nSurf} '\verticesInfo.mat'],'validCellsFinal','newVerticesNeighs2D')
% 
%         neighsCells = cell(1,max(newVerticesNeighs2D(:)));
%         for nCell = 1 : max(newVerticesNeighs2D(:))
%         	 [nRow,nCol]=find(ismember(newVerticesNeighs2D,nCell));
%              cellsNeigh = unique(newVerticesNeighs2D(nRow,:));
%              cellsNeigh = cellsNeigh(cellsNeigh~=nCell);
%              neighsCells{nCell} = cellsNeigh;             
%         end
%         sidesCells = cellfun(@(x) length(x), neighsCells);
%         [polyDist]=calculate_polygon_distribution(sidesCells,validCellsFinal);
%         
%         volumePerim = regionprops3(img3d,'Volume');
%         areaCells = cat(1,volumePerim.Volume);
%         areaValidCells = areaCells(validCellsFinal);
%         
%         if nSurf == 1
%             polyDistApical{nFile} = polyDist(2,:);
%             logNormAreaApical{nFile} = log10(areaValidCells./(mean(areaValidCells)));
%             normAreaApical{nFile} = areaValidCells./(mean(areaValidCells));
%             totalSidesCellsApical{nFile} = sidesCells(validCellsFinal);
%         else
%             polyDistBasal{nFile} = polyDist(2,:);
%             logNormAreaBasal{nFile} = log10(areaValidCells./(mean(areaValidCells)));
%             normAreaBasal{nFile} = areaValidCells./(mean(areaValidCells));
%             totalSidesCellsBasal{nFile} = sidesCells(validCellsFinal);
%         end
%     end
% end
% 
% polyDistBasal = cell2mat(vertcat(polyDistBasal{:}));
% meanPolyDistBasal = mean(polyDistBasal);
% stdPolyDistBasal = std(polyDistBasal);
% 
% polyDistApical = cell2mat(vertcat(polyDistApical{:}));
% meanPolyDistApical = mean(polyDistApical);
% stdPolyDistApical = std(polyDistApical);
% 
% dispersionLogNormAreaBasal = vertcat(logNormAreaBasal{:});
% dispersionLogNormAreaApical = vertcat(logNormAreaApical{:});
% 
% dispersionNormAreaBasal = vertcat(normAreaBasal{:});
% dispersionNormAreaApical = vertcat(normAreaApical{:});
% 
% listPolygons = 3:23;
% 
% relationNormArea_numSidesBasal = [horzcat(totalSidesCellsBasal{:})',dispersionNormAreaBasal];
% relationNormArea_numSidesApical = [horzcat(totalSidesCellsApical{:})',dispersionNormAreaApical];
% 
% uniqSidesBasal = unique(horzcat(totalSidesCellsBasal{:}));
% uniqSidesApical = unique(horzcat(totalSidesCellsApical{:}));
% 
% lewisBasal_NormArea = [uniqSidesBasal;arrayfun(@(x) mean(relationNormArea_numSidesBasal(ismember(relationNormArea_numSidesBasal(:,1),x),2)),uniqSidesBasal);
%     arrayfun(@(x) std(relationNormArea_numSidesBasal(ismember(relationNormArea_numSidesBasal(:,1),x),2)),uniqSidesBasal)];
% 
% lewisApical_NormArea = [uniqSidesApical;arrayfun(@(x) mean(relationNormArea_numSidesApical(ismember(relationNormArea_numSidesApical(:,1),x),2)),uniqSidesApical);
%     arrayfun(@(x) std(relationNormArea_numSidesApical(ismember(relationNormArea_numSidesApical(:,1),x),2)),uniqSidesApical)];
% 
% save(['docs\figuresMathPaper\lewis2D_averagePolygon_AreasDistribution_' date '.mat'],'meanPolyDistBasal','stdPolyDistBasal','meanPolyDistApical','stdPolyDistApical','dispersionLogNormAreaBasal','dispersionLogNormAreaApical','dispersionNormAreaBasal','dispersionNormAreaApical','listPolygons','lewisBasal_NormArea','lewisApical_NormArea')



%%Represent figures
folderSalGland = 'docs\figuresMathPaper\lewis2D_averagePolygon_AreasDistribution_06-Mar-2019.mat';
load(folderSalGland)

%%  figure 1C area distribution
folderTube = 'D:\Pedro\Epithelia3D\InSilicoModels\TubularModel\data\tubularCVT\Data\512x4096_200seeds\polygonDistribution_diag_9.mat';
path2save = 'docs\figuresMathPaper\';
load(folderTube)

h = figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
histogram(dispersionLogNormArea,'BinWidth',0.02,'Normalization','probability','FaceColor',[0,0,0.8],'FaceAlpha',1); 
xlim([-0.8,0.8])
ylim([0,0.1])
title('Voronoi 9 - tube')
xlabel('log10(normalized area)')
ylabel('proportion')
set(gca,'FontSize', 24,'FontName','Helvetica','YGrid','on','Box','off');
print(h,[path2save 'fig1C_voronoi9_' date],'-dtiff','-r300')

h = figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
histogram(dispersionLogNormAreaBasal,'BinWidth',0.02,'Normalization','probability','FaceColor',[0,0.4,0.2],'FaceAlpha',1); 
xlim([-0.8,0.8])
ylim([0,0.1])
title('Basal - Glands')
xlabel('log10(normalized area)')
ylabel('proportion')
set(gca,'FontSize', 24,'FontName','Helvetica','YGrid','on','Box','off');
print(h,[path2save 'fig1C_basalGland_' date],'-dtiff','-r300')

h = figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
histogram(dispersionLogNormAreaApical,'BinWidth',0.02,'Normalization','probability','FaceColor',[0.2,0.8,0],'FaceAlpha',1); 
xlim([-0.8,0.8])
ylim([0,0.1])
title('Apical - Glands')
xlabel('log10(normalized area)')
ylabel('proportion')
set(gca,'FontSize', 24,'FontName','Helvetica','YGrid','on','Box','off');
print(h,[path2save 'fig1C_apicalGland_' date],'-dtiff','-r300')
    
%%  figure 2 Lewis apical basal - glands
% lewisBasal_logNormArea
% lewisApical_logNormArea