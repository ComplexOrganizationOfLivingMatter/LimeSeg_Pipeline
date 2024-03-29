clear all
close all

addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

% %% Extracting and organizing information from segmented 3D Salivary Glands

% files = dir('**/data/Salivary gland_ExtractedVertices_Correct/**/Results/3d_layers_info.mat');
% 
% nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
% nonDiscardedFilesEcadh = cellfun(@(x) contains(lower(x), 'inhibited') == 0, {files.folder});
% 
% files = files((nonDiscardedFiles.*nonDiscardedFilesEcadh)>0);
% 
% surfLayers = {'1','basal'};
% 
% polyDistApical = cell(size(files,1),1);
% polyDistBasal = cell(size(files,1),1);
% logNormAreaApical = cell(size(files,1),1);
% logNormAreaBasal = cell(size(files,1),1);
% normAreaApical = cell(size(files,1),1);
% normAreaBasal = cell(size(files,1),1);
% normVolumeCells = cell(size(files,1),1);
% totalSidesCellsApical = cell(size(files,1),1);
% totalSidesCellsBasal = cell(size(files,1),1);
% totalSidesCells = cell(size(files,1),1);
% neighExchanges = zeros(size(files,1),1);
% percScutoids = zeros(size(files,1),1); 
% 
% 
% for nFile = 1 : size(files,1)
%       
%     for nSurf = 1:2

%         load([files(nFile).folder '\unrolledGlands\gland_SR_' surfLayers{nSurf} '\final3DImg.mat'],'img3d')
%         load([files(nFile).folder '\unrolledGlands\gland_SR_' surfLayers{nSurf} '\verticesInfo.mat'],'validCellsFinal','newVerticesNeighs2D')
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
%         
%         
%         if nSurf == 1
%             polyDistApical{nFile} = polyDist(2,:);
%             logNormAreaApical{nFile} = log10(areaValidCells./(mean(areaValidCells)));
%             normAreaApical{nFile} = areaValidCells./(mean(areaValidCells));
%             totalSidesCellsApical{nFile} = sidesCells(validCellsFinal);
%             neighsCellsApical = neighsCells;
%         else
%             polyDistBasal{nFile} = polyDist(2,:);
%             logNormAreaBasal{nFile} = log10(areaValidCells./(mean(areaValidCells)));
%             normAreaBasal{nFile} = areaValidCells./(mean(areaValidCells));
%             totalSidesCellsBasal{nFile} = sidesCells(validCellsFinal);
%             neighsCellsBasal = neighsCells;
%             
%             neighChanges = cellfun(@(x,y) length(setxor(x,y)),neighsCellsBasal(validCellsFinal),neighsCellsApical(validCellsFinal));
%             numScutoids = cellfun(@(x,y) ~isempty(setxor(x,y)),neighsCellsBasal(validCellsFinal),neighsCellsApical(validCellsFinal));
%             
%             totalSidesCells{nFile} = cellfun(@(x,y) length(unique([x;y])),neighsCellsBasal(validCellsFinal),neighsCellsApical(validCellsFinal));
% 
%             
%             neighExchanges(nFile) = sum(neighChanges)/length(validCellsFinal);
%             percScutoids(nFile) = sum(numScutoids)/length(validCellsFinal);
%             
%             
%         end
%     end
%     load([files(nFile).folder '\3d_layers_info.mat'],'labelledImage')
%     
%     volumeCells = regionprops3(labelledImage,'Volume');
%     volumeCells = cat(1,volumeCells.Volume);
%     volumeValidCells = volumeCells(validCellsFinal);
%     normVolumeCells{nFile} = volumeValidCells./(mean(volumeValidCells));
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
% dispersionNormVolume = vertcat(normVolumeCells{:});
% 
% listPolygons = 3:23;
% 
% relationNormArea_numSidesBasal = [horzcat(totalSidesCellsBasal{:})',dispersionNormAreaBasal];
% relationNormArea_numSidesApical = [horzcat(totalSidesCellsApical{:})',dispersionNormAreaApical];
% relationNormVol_numSidesVolume = [horzcat(totalSidesCells{:})',dispersionNormVolume];
% 
% 
% %relationNsidesApical-basal
% totalSidesCells = horzcat(totalSidesCells{:})';
% totalSidesTotalApical = arrayfun(@(x) totalSidesCells(horzcat(totalSidesCellsApical{:})' == x),4:9,'UniformOutput',false);
% totalSidesTotalBasal = arrayfun(@(x) totalSidesCells(horzcat(totalSidesCellsBasal{:})' == x),4:9,'UniformOutput',false);
% 
 
% uniqSidesBasal = unique(horzcat(totalSidesCellsBasal{:}));
% uniqSidesApical = unique(horzcat(totalSidesCellsApical{:}));
% uniqSides3D = unique(totalSidesCells);
% 
% 
% lewisBasal_NormArea = [uniqSidesBasal;arrayfun(@(x) mean(relationNormArea_numSidesBasal(ismember(relationNormArea_numSidesBasal(:,1),x),2)),uniqSidesBasal);
%     arrayfun(@(x) std(relationNormArea_numSidesBasal(ismember(relationNormArea_numSidesBasal(:,1),x),2)),uniqSidesBasal)];
% 
% lewisApical_NormArea = [uniqSidesApical;arrayfun(@(x) mean(relationNormArea_numSidesApical(ismember(relationNormArea_numSidesApical(:,1),x),2)),uniqSidesApical);
%     arrayfun(@(x) std(relationNormArea_numSidesApical(ismember(relationNormArea_numSidesApical(:,1),x),2)),uniqSidesApical)];
% 
% lewis3D_NormVol = [uniqSides3D';arrayfun(@(x) mean(relationNormVol_numSidesVolume(ismember(relationNormVol_numSidesVolume(:,1),x),2)),uniqSides3D)';
%     arrayfun(@(x) std(relationNormVol_numSidesVolume(ismember(relationNormVol_numSidesVolume(:,1),x),2)),uniqSides3D)'];
% 
% 
% meanScutoids = mean(percScutoids);
% stdScutoids = std(percScutoids);
% meanNeighExchanges = mean(neighExchanges);
% stdNeighExchanges = std(neighExchanges);
% 
% mkdir('docs\figuresMathPaper\');

% save(['docs\figuresMathPaper\lewis2D_3D_averagePolygon_AreasDistribution_' date '.mat'],'meanScutoids','stdScutoids','meanNeighExchanges','stdNeighExchanges','meanPolyDistBasal','stdPolyDistBasal','meanPolyDistApical','stdPolyDistApical','dispersionLogNormAreaBasal','dispersionLogNormAreaApical','dispersionNormAreaBasal','dispersionNormAreaApical','dispersionNormVolume','listPolygons','lewisBasal_NormArea','lewisApical_NormArea','lewis3D_NormVol','totalSidesTotalApical','totalSidesTotalBasal')



%% Represent paper figures (LEWIS 2D and 3D, and 'poor get richer')

% folderSalGland = 'docs\figuresMathPaper\lewis2D_3D_averagePolygon_AreasDistribution_29-Mar-2019.mat';

folderSalGland = ['docs\figuresMathPaper\lewis2D_3D_averagePolygon_AreasDistribution_09-Nov-2019.mat'];

load(folderSalGland)
lewis3D_glands = lewis3D_NormVol;
lewisApical_glands = lewisApical_NormArea;
lewisBasal_glands = lewisBasal_NormArea;

%%  figure area distribution
VorN = 8;

if VorN == 8
    colorPlotTubApi = [79,209,255]./255;
    colorPlotT_175 = [0,112/255,192/255];
    colorPlotT_4 = [27/255,39/255,201/255];
else
    colorPlotTubApi = [210,210,210]./255;
    colorPlotT_175 = [190,190,190]./255;
    colorPlotT_4 = [170,170,170]./255;
end

colorPlotApiG = [0.2,0.8,0];
colorPlotBasG = [0,0.4,0.2];
colorPlotGland = round(mean([[0.2,0.8,0];[0,0.4,0.2]]),2);


%% load tubes data
nRealizations = 20;
H = 512;
W = 4096;
nSeeds = 200;

folderTube = ['..\Epithelia3D\InSilicoModels\TubularModel\data\tubularCVT\Data\' num2str(H) 'x' num2str(W) '_' num2str(nSeeds) 'seeds\'];
if ~exist([folderTube 'polygonDistribution_diag_' num2str(VorN) '.mat'],'file')
    calculatePolDist(H,W,nSeeds,VorN,nRealizations)
end
load([folderTube 'polygonDistribution_diag_' num2str(VorN) '.mat'])
lewisApical_tube = lewis_NormArea; 
path2save = 'docs\figuresMathPaper\';



%% tube basal - 1.75
folderTube = strrep(folderTube,'tubularCVT\Data','tubularVoronoiModel\expansion');
folderTube = [folderTube 'diagram' num2str(VorN) '_Markov\'];
load([folderTube 'polygonDistribution_diag_' num2str(VorN) 'sr1.75_volume.mat']);
apicalSidesCellsTube = apicalSidesCells;
basalSidesCellsTube175 = basalSidesCells;
totalSidesCellsTube175 = totalSidesCells;
neighsApical=neighsApicalRealizations;
neighsBasal175=neighsBasalRealizations;
validCells_tube_175 = validCellsRealizations;

lewisBasal_tube_175 = lewis_NormArea; 
lewis3D_tube_175 = lewis3D_volNorm;
validCells_tube_175 = validCellsRealizations;

%% tube basal - 4
load([folderTube 'polygonDistribution_diag_' num2str(VorN) 'sr4_volume.mat'])

basalSidesCellsTube4 = basalSidesCells;
totalSidesCellsTube4 = totalSidesCells;
neighsBasal4=neighsBasalRealizations;
validCells_tube_4 = validCellsRealizations;

lewisBasal_tube_4 = lewis_NormArea; 
lewis3D_tube_4 = lewis3D_volNorm;




%%Group and apply valid cells
validCells_global = cellfun(@(x,y) intersect(x,y),validCells_tube_175,validCells_tube_4,'UniformOutput',false);

for nRea = 1:nRealizations
    
    neighsAp = neighsApical{nRea};
    neighsApical{nRea}=neighsAp(validCells_global{nRea});
    neighs175 = neighsBasal175{nRea};
    neighsBasal175{nRea}=neighs175(validCells_global{nRea});
    neighs4 = neighsBasal4{nRea};
    neighsBasal4{nRea}=neighs4(validCells_global{nRea});
    
    sidesAp = apicalSidesCellsTube{nRea};
    apicalSidesCellsTube{nRea} = sidesAp(validCells_global{nRea});
    sidesBas175 = basalSidesCellsTube175{nRea};
    basalSidesCellsTube175{nRea} = sidesBas175(validCells_global{nRea});
    sidesBas4 = basalSidesCellsTube4{nRea};
    basalSidesCellsTube4{nRea} = sidesBas4(validCells_global{nRea});
    
    sidestot175 = totalSidesCellsTube175{nRea};
    totalSidesCellsTube175{nRea} = sidestot175(validCells_global{nRea});
    sidestot4 = totalSidesCellsTube4{nRea};
    totalSidesCellsTube4{nRea} = sidestot4(validCells_global{nRea});
    
    totalSidesCellsTubeFrom175toApical{nRea} = cellfun(@(x,y) length(unique([x;y])),neighsBasal175{nRea},neighsApical{nRea});
    totalSidesCellsTubeFrom4toApical{nRea} = cellfun(@(x,y,z) length(unique([x;y;z])),neighsBasal4{nRea},neighsBasal175{nRea},neighsApical{nRea});
end


% %% Fig LEWIS 2D
%     h = figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
%     markerSiz = 20;
%     linWid = 2;
%     
%     indValApiG = ismember(lewisApical_glands(1,:),5:7);
%     pA = polyfit(lewisApical_glands(1,indValApiG), lewisApical_glands(2,indValApiG),1);
%     f = polyval(pA,[min(lewisApical_glands(1,indValApiG))-0.5, lewisApical_glands(1,indValApiG), max(lewisApical_glands(1,indValApiG))+0.5]); 
%     hold on
%     plot([min(lewisApical_glands(1,indValApiG))-0.5, lewisApical_glands(1,indValApiG), max(lewisApical_glands(1,indValApiG))+0.5],f,'Color',colorPlotApiG,'LineWidth',linWid)
%     mdlApical_glands = polyfitn(lewisApical_glands(1,indValApiG), lewisApical_glands(2,indValApiG),1);
%    
%     indValBasG = ismember(lewisBasal_glands(1,:),5:7);
%     pB = polyfit(lewisBasal_glands(1,indValBasG), lewisBasal_glands(2,indValBasG),1);
%     f = polyval(pB,[min(lewisBasal_glands(1,indValBasG))-0.5, lewisBasal_glands(1,indValBasG), max(lewisBasal_glands(1,indValBasG))+0.5]); 
%     hold on
%     plot([min(lewisBasal_glands(1,indValBasG))-0.5, lewisBasal_glands(1,indValBasG), max(lewisBasal_glands(1,indValBasG))+0.5],f,'Color',colorPlotBasG,'LineWidth',linWid)
%     mdlBasal_glands = polyfitn(lewisBasal_glands(1,indValBasG), lewisBasal_glands(2,indValBasG),1);
%    
%     indValTubApi = ismember(lewisApical_tube(1,:),4:8);
%     pTa = polyfit(lewisApical_tube(1,indValTubApi), lewisApical_tube(2,indValTubApi),1);
%     f = polyval(pTa,[min(lewisApical_tube(1,indValTubApi))-0.5, lewisApical_tube(1,indValTubApi), max(lewisApical_tube(1,indValTubApi))+0.5]); 
%     plot([min(lewisApical_tube(1,indValTubApi))-0.5, lewisApical_tube(1,indValTubApi), max(lewisApical_tube(1,indValTubApi))+0.5],f,'Color',colorPlotTubApi,'LineWidth',linWid)
%     mdlT_apical = polyfitn(lewisApical_tube(1,indValTubApi), lewisApical_tube(2,indValTubApi),1);
% 
%     indValTub_175 = ismember(lewisBasal_tube_175(1,:),4:8);
%     pTb175 = polyfit(lewisBasal_tube_175(1,indValTub_175), lewisBasal_tube_175(2,indValTub_175),1);
%     f = polyval(pTb175,[min(lewisBasal_tube_175(1,indValTub_175))-0.5, lewisBasal_tube_175(1,indValTub_175), max(lewisBasal_tube_175(1,indValTub_175))+0.5]); 
%     plot([min(lewisBasal_tube_175(1,indValTub_175))-0.5, lewisBasal_tube_175(1,indValTub_175), max(lewisBasal_tube_175(1,indValTub_175))+0.5],f,'Color',colorPlotT_175,'LineWidth',linWid)
%     mdlT_175 = polyfitn(lewisBasal_tube_175(1,indValTub_175), lewisBasal_tube_175(2,indValTub_175),1);
%   
%     indValTub_4 = ismember(lewisBasal_tube_4(1,:),4:8);
%     pTb4 = polyfit(lewisBasal_tube_4(1,indValTub_4), lewisBasal_tube_4(2,indValTub_4),1);
%     f = polyval(pTb4,[min(lewisBasal_tube_4(1,indValTub_4))-0.5, lewisBasal_tube_4(1,indValTub_4), max(lewisBasal_tube_4(1,indValTub_4))+0.5]); 
%     plot([min(lewisBasal_tube_4(1,indValTub_4))-0.5, lewisBasal_tube_4(1,indValTub_4), max(lewisBasal_tube_4(1,indValTub_4))+0.5],f,'Color',colorPlotT_4,'LineWidth',linWid)
%     mdlT_4 = polyfitn(lewisBasal_tube_4(1,indValTub_4), lewisBasal_tube_4(2,indValTub_4),1);
%     
%     plot(lewisApical_glands(1,indValApiG),lewisApical_glands(2,indValApiG),'.','Color',colorPlotApiG,'MarkerSize',markerSiz)
%     plot(lewisBasal_glands(1,indValBasG),lewisBasal_glands(2,indValBasG),'.','Color',colorPlotBasG,'MarkerSize',markerSiz)
%     plot(lewisApical_tube(1,indValTubApi),lewisApical_tube(2,indValTubApi),'.','Color',colorPlotTubApi,'MarkerSize',markerSiz)
%     plot(lewisBasal_tube_175(1,indValTub_175),lewisBasal_tube_175(2,indValTub_175),'.','Color',colorPlotT_175,'MarkerSize',markerSiz)
%     plot(lewisBasal_tube_4(1,indValTub_4),lewisBasal_tube_4(2,indValTub_4),'.','Color',colorPlotT_4,'MarkerSize',markerSiz)
% 
%     title('Lewis 2D')
%     ylabel('Area normalized')
%     xlabel('neighbours total')
%     xlim([3 9])
%     ylim([0.3 1.6])
%     yticks([0.5 0.75 1 1.25 1.5])
%     set(gca,'FontSize', 24,'FontName','Helvetica','YGrid','on','TickDir','out','Box','off');
%     
%     print(h,[path2save 'fig_Lewis2D_Vor' num2str(VorN) '_' date '_noLegend'],'-dtiff','-r300')
%     
%     legend({['fit R2 ' num2str(round(mdlApical_glands.R2,4)) ' - RMSE '  num2str(round(mdlApical_glands.RMSE,4)) ' - slope ' num2str(round(pA(1),3))],...
%             ['fit R2 ' num2str(round(mdlBasal_glands.R2,4)) ' - RMSE '  num2str(round(mdlBasal_glands.RMSE,4)) ' - slope ' num2str(round(pB(1),3))],...
%             ['fit R2 ' num2str(round(mdlT_apical.R2,4)) ' - RMSE '  num2str(round(mdlT_apical.RMSE,4)) ' - slope ' num2str(round(pTa(1),3))],...
%             ['fit R2 ' num2str(round(mdlT_175.R2,4)) ' - RMSE '  num2str(round(mdlT_175.RMSE,4)) ' - slope ' num2str(round(pTb175(1),3))],...
%             ['fit R2 ' num2str(round(mdlT_4.R2,4)) ' - RMSE '  num2str(round(mdlT_4.RMSE,4)) ' - slope ' num2str(round(pTb4(1),3))],'Glands apical',...
%             'Glands basal',['Voronoi tube ' num2str(VorN) ' apical'],['Voronoi tube ' num2str(VorN) ' basal 1.8'],['Voronoi tube ' num2str(VorN) ' basal 4']},'FontSize',12)
%     hold off
%     print(h,[path2save 'fig_Lewis2D_Vor_' num2str(VorN) '_' date],'-dtiff','-r300')
%     savefig(h,[path2save 'fig_Lewis2D_Vor_' num2str(VorN) '_' date])


% %% Fig LEWIS 3D
%     h = figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
%     markerSiz = 20;
%     linWid = 2;   
%     hold on
% 
%     indValT175 = ismember(lewis3D_tube_175(1,:),5:8);
%     pT175 = polyfit(lewis3D_tube_175(1,indValT175), lewis3D_tube_175(2,indValT175),1);
%     f = polyval(pT175,[min(lewis3D_tube_175(1,indValT175))-0.5, lewis3D_tube_175(1,indValT175), max(lewis3D_tube_175(1,indValT175))+0.5]); 
%     plot([min(lewis3D_tube_175(1,indValT175))-0.5, lewis3D_tube_175(1,indValT175), max(lewis3D_tube_175(1,indValT175))+0.5],f,'Color',colorPlotT_175,'LineWidth',linWid)
%     mdlT_175 = polyfitn(lewis3D_tube_175(1,indValT175), lewis3D_tube_175(2,indValT175),1);
%     
%     indValT4 = ismember(lewis3D_tube_4(1,:),7:11);
%     pT4 = polyfit(lewis3D_tube_4(1,indValT4), lewis3D_tube_4(2,indValT4),1);
%     f = polyval(pT4,[min(lewis3D_tube_4(1,indValT4))-0.5, lewis3D_tube_4(1,indValT4), max(lewis3D_tube_4(1,indValT4))+0.5]); 
%     plot([min(lewis3D_tube_4(1,indValT4))-0.5, lewis3D_tube_4(1,indValT4), max(lewis3D_tube_4(1,indValT4))+0.5],f,'Color',colorPlotT_4,'LineWidth',linWid)
%     mdlT_4 = polyfitn(lewis3D_tube_4(1,indValT4), lewis3D_tube_4(2,indValT4),1);
% 
%     indValGl = ismember(lewis3D_glands(1,:),6:8);
%     pTSG = polyfit(lewis3D_glands(1,indValGl), lewis3D_glands(2,indValGl),1);
%     f = polyval(pTSG,[min(lewis3D_glands(1,indValGl))-0.5, lewis3D_glands(1,indValGl), max(lewis3D_glands(1,indValGl))+0.5]); 
%     plot([min(lewis3D_glands(1,indValGl))-0.5, lewis3D_glands(1,indValGl), max(lewis3D_glands(1,indValGl))+0.5],f,'Color',colorPlotGland,'LineWidth',linWid)
%     mdlGlands = polyfitn(lewis3D_glands(1,indValGl), lewis3D_glands(2,indValGl),1);
%     
%     plot(lewis3D_tube_175(1,indValT175),lewis3D_tube_175(2,indValT175),'.','Color',colorPlotT_175,'MarkerSize',markerSiz)
%     plot(lewis3D_tube_4(1,indValT4),lewis3D_tube_4(2,indValT4),'.','Color',colorPlotT_4,'MarkerSize',markerSiz)    
%     plot(lewis3D_glands(1,indValGl),lewis3D_glands(2,indValGl),'.','Color',colorPlotGland,'MarkerSize',markerSiz)
% 
% 
%     title('Lewis 3D')
%     ylabel('Volume normalized')
%     xlabel('neighbours total')
%     xlim([4 12])
%     ylim([0.3 1.5])
%     yticks([0.5 0.75 1 1.25 1.5])
%     set(gca,'FontSize', 24,'FontName','Helvetica','YGrid','on','TickDir','out','Box','off');
%     
% %     preD = predint(mdlT8_4,unique([lewis3D_tube8_175(1,indValT8175),lewis3D_tube8_4(1,indValT84)]),0.95,'observation','off');
% %     plot([surfRatios max(surfRatios)+1],preD,'--','Color',colorPlot)
% %     
%     print(h,[path2save 'fig_Lewis3D_Vor_' num2str(VorN) '_' date '_noLegend'],'-dtiff','-r300')
% 
%     legend({['fit R2 ' num2str(round(mdlT_175.R2,4)) ' - RMSE '  num2str(round(mdlT_175.RMSE,4)) ' - slope ' num2str(round(pT175(1),3))],...
%         ['fit R2 ' num2str(round(mdlT_4.R2,4)) ' - RMSE '  num2str(round(mdlT_4.RMSE,4)) ' - slope ' num2str(round(pT4(1),3))],...
%         ['fit R2 ' num2str(round(mdlGlands.R2,4)) ' - RMSE '  num2str(round(mdlGlands.RMSE,4)) ' - slope ' num2str(round(pTSG(1),3))]...
%         ['Voronoi tube ' num2str(VorN) ' - 1.8'],['Voronoi tube ' num2str(VorN) ' - 4'],'Glands'},'FontSize',12)
%     hold off
%     
%     print(h,[path2save 'fig_Lewis3D_Vor_' num2str(VorN) '_' date],'-dtiff','-r300')
%     savefig(h,[path2save 'fig_Lewis3D_Vor_' num2str(VorN) '_' date])


%%  figure Relation apical - basal nSides. 'Poor get richer'

%from apical to basal
sidesTotalApical = arrayfun(@(x) totalSidesTotalApical{x-3}-x,4:9,'UniformOutput',false);
totalSidesCellsTube175 = horzcat(totalSidesCellsTube175{:});

totalSidesCellsTube4 = cellfun(@(x,y) x(y), totalSidesCellsTube4,validCells_tube_4,'UniformOutput',false);
apicalSidesCellsTube4 = cellfun(@(x,y) x(y), apicalSidesCellsTube4,validCells_tube_4,'UniformOutput',false);
totalSidesCellsTube4 = horzcat(totalSidesCellsTube4{:});

sidesTotalTub175 = arrayfun(@(x) totalSidesCellsTube175(horzcat(apicalSidesCellsTube{:})' == x)-x,4:9,'UniformOutput',false);
sidesTotalTub4 = arrayfun(@(x) totalSidesCellsTube4(horzcat(apicalSidesCellsTube{:})' == x)-x,4:9,'UniformOutput',false);


%from basal to apical
sidesTotalBasal = arrayfun(@(x) totalSidesTotalBasal{x-3}-x,4:9,'UniformOutput',false);
totalSidesCellsFrom175toApical = horzcat(totalSidesCellsTubeFrom175toApical{:});
totalSidesCellsFrom4toApical = horzcat(totalSidesCellsTubeFrom4toApical{:});

sidesTotalTub175toApical = arrayfun(@(x) totalSidesCellsFrom175toApical(horzcat(basalSidesCellsTube175{:})' == x)-x,4:9,'UniformOutput',false);
sidesTotalTub4toApical = arrayfun(@(x) totalSidesCellsFrom4toApical(horzcat(basalSidesCellsTube4{:})' == x)-x,4:9,'UniformOutput',false);

%selectMethod {from apical to basal -> 0, from basal to apical ->1}
selectMethod=1;
if selectMethod == 1
    sidesTotalGenericInitGl = sidesTotalBasal;
    sidesTotalGenericTubeStep1 = sidesTotalTub175toApical;
    sidesTotalGenericTubeStep2 = sidesTotalTub4toApical;
else
    sidesTotalGenericInitGl = sidesTotalApical;
    sidesTotalGenericTubeStep1 = sidesTotalTub175;
    sidesTotalGenericTubeStep2 = sidesTotalTub4;
end
h = figure('units','normalized','outerposition',[0 0 1 1],'Visible','on');
hold on

% nInitGland = cellfun(@(x,y) ones(size(x))*y,sidesTotalGenericInitGl,[{4},{5},{6},{7},{8},{9}],'UniformOutput',false);
% relGland = [vertcat(nInitGland{:}),vertcat(sidesTotalGenericInitGl{:})];
% 
% nInitTube1 = cellfun(@(x,y) (ones(size(x))*y),sidesTotalGenericTubeStep1,[{4},{5},{6},{7},{8},{9}],'UniformOutput',false);
% relTube1 = [horzcat(nInitTube1{:});horzcat(sidesTotalGenericTubeStep1{:})]';
% 
% nInitTube2 = cellfun(@(x,y) (ones(size(x))*y),sidesTotalGenericTubeStep2,[{4},{5},{6},{7},{8},{9}],'UniformOutput',false);
% relTube2 = [horzcat(nInitTube2{:});horzcat(sidesTotalGenericTubeStep2{:})]';

x = [sidesTotalGenericInitGl{1};sidesTotalGenericTubeStep1{1}';sidesTotalGenericTubeStep2{1}';...
    sidesTotalGenericInitGl{2};sidesTotalGenericTubeStep1{2}';sidesTotalGenericTubeStep2{2}';...
    sidesTotalGenericInitGl{3};sidesTotalGenericTubeStep1{3}';sidesTotalGenericTubeStep2{3}';...
    sidesTotalGenericInitGl{4};sidesTotalGenericTubeStep1{4}';sidesTotalGenericTubeStep2{4}';...
    sidesTotalGenericInitGl{5};sidesTotalGenericTubeStep1{5}';sidesTotalGenericTubeStep2{5}'];

g = [zeros(length(sidesTotalGenericInitGl{1}), 1); ones(length(sidesTotalGenericTubeStep1{1}), 1);2*ones(length(sidesTotalGenericTubeStep2{1}), 1);...
    3*ones(length(sidesTotalGenericInitGl{2}), 1); 4*ones(length(sidesTotalGenericTubeStep1{2}), 1);5*ones(length(sidesTotalGenericTubeStep2{2}), 1);...
    6*ones(length(sidesTotalGenericInitGl{3}), 1); 7*ones(length(sidesTotalGenericTubeStep1{3}), 1);8*ones(length(sidesTotalGenericTubeStep2{3}), 1);...
    9*ones(length(sidesTotalGenericInitGl{4}), 1); 10*ones(length(sidesTotalGenericTubeStep1{4}), 1);11*ones(length(sidesTotalGenericTubeStep2{4}), 1);...
    12*ones(length(sidesTotalGenericInitGl{5}), 1); 13*ones(length(sidesTotalGenericTubeStep1{5}), 1);14*ones(length(sidesTotalGenericTubeStep2{5}), 1)];

%G = iosr.statistics.boxPlot(x',g');

H = notBoxPlot(x,g,'markMedian',true,'jitter', 0.3);

set([H.data],'MarkerSize',4,...
    'markerFaceColor','none',...
    'markerEdgeColor', 'none')
%mean line color
set([H.mu],'color','k')

for ii=1:length(H)
    set(H(ii).perc1,'FaceColor','none',...
                   'EdgeColor','k','lineStyle','-')
          
    set(H(ii).sd,'FaceColor','none',...
                   'EdgeColor',[0.6 0.6 0.6],'lineStyle',':')        
    set(H(ii).sd1,'Color',[0.6 0.6 0.6])     
    set(H(ii).sd2,'Color',[0.6 0.6 0.6])            
%     set(H(ii).sdPtch,'FaceColor','none',...
%                    'EdgeColor','k')
%     set(H(ii).semPtch,'FaceColor','none',...
%                    'EdgeColor','k')
end

typeOfPoly = 4:8;
for nSideAp = 1:10 %(4:8)
    switch nSideAp

        case {1,4,7,10,13}
            freqSides = sidesTotalGenericInitGl{ceil(nSideAp/3)};
            c = [0.2,0.8,0];
        case {2,5,8,11,14}
            freqSides = sidesTotalGenericTubeStep1{ceil(nSideAp/3)};
            c = [0,112/255,192/255];
        otherwise
            freqSides = sidesTotalGenericTubeStep2{ceil(nSideAp/3)};
            c = [27/255,39/255,201/255];
    end
    
    typeCell = typeOfPoly(ceil(nSideAp/2));
    
    numbers=unique(freqSides);       %list of elements
    countN=arrayfun(@(x) sum(ismember(freqSides,x)),numbers); 
    countP = countN./sum(countN);
    try
        for nUniqSides = 1:length(numbers)
            hold on
%             scatter(nSideAp-1, numbers(nUniqSides),4000*countP(nUniqSides), 'filled','MarkerFaceColor',c, 'MarkerEdgeColor', 'k','MarkerFaceAlpha',0.5)
            text(nSideAp-1,numbers(nUniqSides),num2str(countN(nUniqSides)),'HorizontalAlignment','center','FontSize',12)
        end
    catch
        
    end
end


% xlim([-1 30])
% xticks([-1:30])
ylim([-1 6])
title('relation apical sides - added neigh')
xlabel('number of sides - apical')
ylabel('number gain total')
set(gca,'FontSize', 24,'FontName','Helvetica','YGrid','on','TickDir','out','Box','off');

xticklabels({'4','5','6','7','8'})

savefig(h,[path2save 'fig_PoorGetRicherFromBasal_' date]);
print(h,[path2save 'fig_PoorGetRicherFromBasal_' date],'-dtiff','-r300')










