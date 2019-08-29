%% Messy file to create the figures of the article, obtaining several measurements which were not captured directly.
% Use it at your own risk.
clear all
close all

addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files_WT = dir('**/data/Salivary gland_ExtractedVertices_Correct/**/Results/3d_layers_info.mat');
files_Ecadhi = dir('**/data/Salivary gland/E-cadh Inhibited/**/Results/3d_layers_info.mat');

files = files_Ecadhi;
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
files = files(nonDiscardedFiles);

%resultsFileName = '3d_layers_info.mat';
%resultsFileName = 'glandDividedInSurfaceRatios.mat';
%resultsFileName = 'glandDividedInSurfaceRatios_PredefinedSR.mat';
resultsFileName = 'glandDividedInSurfaceRatios_AllUnrollFeatures.mat';

minNumberOfSurfaceRatios = 5; %7WT; 5Ecadhi
%namesSR = arrayfun(@(x) ['sr' strrep(num2str(x),'.','_')],1:numberOfSurfaceRatios,'UniformOutput', false);
%namesSR = {'sr1' 'sr2'};

steps = 2.5/(7-1);
surfaceRatiosExtrapolatedFrom3D = 1:steps:((steps*(7-1))+1);

totalDataBasal = [];
totalDataAccum = [];

realSRBasal = [];
realSRAccum = [];
meanVolumeMicronsPerGland = zeros(length(files),2);



% if ~exist('Results/salivaryGland_Info_20_05_2019.mat','file')
    for numFile = 1:length(files)
        files(numFile).folder

        load(fullfile(files(numFile).folder, files(numFile).name));
        load(fullfile(files(numFile).folder, 'valid_cells.mat'));
        if exist(fullfile(files(numFile).folder, resultsFileName), 'file')
            load(fullfile(files(numFile).folder, resultsFileName))
        else
            continue
        end
        load([files(numFile).folder '\unrolledGlands\gland_SR_basal\final3DImg.mat'],'img3dComplete')
        
        %% Calculate final variables
        %Volume
        volume3d = regionprops3(img3dComplete,'Volume');
        volume3d = cat(1,volume3d.Volume);
        meanVolumeMicronsPerGland(numFile,1) = mean(volume3d(validCells)); 
        meanVolumeMicronsPerGland(numFile,2) = std(volume3d(validCells));
        
        %Surface ratio: final 3D, final 2D, selectedSR 3D, selectedSR 2D
        meanSR(numFile, 1) = infoPerSurfaceRatio{end, 2};
        meanSR(numFile, 2) = infoPerSurfaceRatio{end, 7};
        meanSR(numFile, 3) = infoPerSurfaceRatio{minNumberOfSurfaceRatios, 2};
        meanSR(numFile, 4) = infoPerSurfaceRatio{minNumberOfSurfaceRatios, 7};
        
        %% Calculate variables per surface
        %Initialize all variables
        numberOfSurfaceRatios = size(infoPerSurfaceRatio, 1);
        neighsSurface = cell(numberOfSurfaceRatios,1);
        neighsAccumSurfaces = cell(numberOfSurfaceRatios,1);
        percentageScutoids = cell(numberOfSurfaceRatios, 1);
        apicoBasalTransitions = cell(numberOfSurfaceRatios, 1);
        areaCells = cell(numberOfSurfaceRatios,1);
        volumes = cell(numberOfSurfaceRatios,1);
        
        infoOfCells = infoPerSurfaceRatio{1, 4};
        infoOfCells = infoOfCells{:};
        
        %Apical sides
        neighsSurface{1} = neighboursOfAllSurfaces{1};
        neighsAccumSurfaces{1} = neighboursOfAllSurfaces{1};
        percentageScutoids{1} = cellfun(@(x, y) ~isequal(x,y), neighsSurface{1}, neighsAccumSurfaces{1});
        numLostNeighsAccum{1} = cell(size(neighboursOfAllSurfaces{1}));
        numWonNeighsAccum{1} = cell(size(neighboursOfAllSurfaces{1}));
        areaCells(1) = {infoOfCells.Basal_area};
        volumes(1) = {infoOfCells.Volume};
        
        for numberSR = 2:numberOfSurfaceRatios
            neighsSurface{numberSR} = neighboursOfAllSurfaces{numberSR};
            neighsAccumSurfaces{numberSR} = cellfun(@(x,y) unique([x;y]),neighsAccumSurfaces{numberSR-1},neighsSurface{numberSR},'UniformOutput',false);
            percentageScutoids{numberSR} = cellfun(@(x, y) ~isempty(setxor(x,y)), neighsSurface{1}, neighsSurface{numberSR});

            lostNeigh = cellfun(@(x, y) setdiff(x,y), neighsAccumSurfaces{numberSR-1}, neighsSurface{numberSR}, 'UniformOutput',false);
            wonNeigh = cellfun(@(x, y) setdiff(y, x), neighsAccumSurfaces{numberSR-1}, neighsAccumSurfaces{numberSR}, 'UniformOutput',false);

            numLostNeighsAccum{numberSR} = cellfun(@(x,y) unique([x;y]),lostNeigh,numLostNeighsAccum{numberSR-1},'UniformOutput',false);
            numWonNeighsAccum{numberSR} = cellfun(@(x,y) unique([x;y]),wonNeigh,numWonNeighsAccum{numberSR-1},'UniformOutput',false);

            apicoBasalTransitions{numberSR} = cellfun(@(x,y) length(([x;y])),numLostNeighsAccum{numberSR},numWonNeighsAccum{numberSR});  

            infoOfCells = infoPerSurfaceRatio{numberSR, 4};
            infoOfCells = infoOfCells{:};
            areaCells(numberSR) = {infoOfCells.Basal_area};
            volumes(numberSR) = {infoOfCells.Volume};
        end

        areaCellsPerSurfaceRealization = cat(2,areaCells{:});
        volumePerSurfaceRealization = cat(2,volumes{:});
        neighsSurface = cat(1,neighsSurface{:})';
        neighsAccumSurfaces = cat(1,neighsAccumSurfaces{:})';
        percentageScutoids = cat(1,percentageScutoids{:})';
        apicoBasalTransitions = cat(1,apicoBasalTransitions{:})';
        numLostNeighsAccum = cat(1, numLostNeighsAccum{:})';
        numWonNeighsAccum = cat(1, numWonNeighsAccum{:})';

        numNeighOfNeighPerSurfacesRealization = getNumNeighsOfNeighs(neighsSurface);
        numNeighOfNeighAccumPerSurfacesRealization = getNumNeighsOfNeighs(neighsAccumSurfaces);

        %%
        meanNumNeighPerSurfaceRealization = mean(numNeighAccumPerSurfacesRealization(validCells, :), 1);
        numCells = repmat(length(validCells), 1, size(meanNumNeighPerSurfaceRealization, 2));
        stdNumNeighPerSurfaceRealization = std(numNeighAccumPerSurfacesRealization(validCells, :), 1);
        totalAreaPerSR = sum(areaCellsPerSurfaceRealization(validCells, :));

        mean_PercScutoids = mean(percentageScutoids(validCells, :), 1);
        std_PercScutoids = std(percentageScutoids(validCells, :), 1);

        mean_apicoBasalTransitions = [0, mean(apicoBasalTransitions(validCells, :), 1)];
        std_apicoBasalTransitions = [0, std(apicoBasalTransitions(validCells, :), 1)];

        mean_apicoBasalTransitionsPerGland{numFile} = mean_apicoBasalTransitions;
        mean_neighsAccum{numFile} = meanNumNeighPerSurfaceRealization;
        mean_Scutoids{numFile} = mean_PercScutoids;
        surfaceRatiosPerFile{numFile} = infoPerSurfaceRatio.SR2D';
    %     
    %     surfaceRatioOfGland_real = vertcat(infoPerSurfaceRatio{:, 7})'; 
    %     totalPartitions = 10;
    %     initialPartitions = (1:(totalPartitions-1))/totalPartitions;
    %     surfaceRatioOfGland = surfaceRatioOfGland_real;
    %     surfaceRatioOfGland(2:10) = initialPartitions * (surfaceRatioOfGland_real(end) - 1) + 1;

        %surfaceRatioOfGland = vertcat(infoPerSurfaceRatio{:, 2})';
        surfaceRatioOfGland = surfaceRatiosExtrapolatedFrom3D;

        mean_PercScutoids_basal(numFile, 1) = mean_PercScutoids(end);
        mean_PercScutoids_basal(numFile, 2) = std_PercScutoids(end);
        mean_PercScutoids_basal(numFile, 3) = infoPerSurfaceRatio{end, 7};
        mean_apicoBasalTransitions_final(numFile, 1) = mean_apicoBasalTransitions(end);
        mean_apicoBasalTransitions_final(numFile, 2) = std_apicoBasalTransitions(end);

        neighbours_apical(numFile, 1) = meanNumNeighPerSurfaceRealization(1);
        neighbours_apical(numFile, 2) = stdNumNeighPerSurfaceRealization(1);
        neighbours_basal(numFile, 1) = mean(numNeighPerSurfaceRealization(validCells, end));
        neighbours_basal(numFile, 2) = std(numNeighPerSurfaceRealization(validCells, end));
        neighbours_total(numFile, 1) = meanNumNeighPerSurfaceRealization(end);
        neighbours_total(numFile, 2) = stdNumNeighPerSurfaceRealization(end);

        %namesSR = surfaceRatioOfGland(1:minNumberOfSurfaceRatios);
        namesSR = surfaceRatiosExtrapolatedFrom3D;
        namesSR = arrayfun(@(x) ['sr' strrep(num2str(x),'.','_')], namesSR, 'UniformOutput', false);

        [polygon_distribution_apical_actual] = calculate_polygon_distribution( numNeighPerSurfaceRealization(:, 1), validCells);
        [polygon_distribution_basal_actual] = calculate_polygon_distribution( numNeighPerSurfaceRealization(:, end), validCells);
        polygon_distribution{numFile} = [polygon_distribution_apical_actual; polygon_distribution_basal_actual(2, :)];

        polygon_distribution_apical(numFile, :) = [polygon_distribution_apical_actual{2, :}];
        polygon_distribution_basal(numFile, :) = [polygon_distribution_basal_actual{2, :}];

        numCells_Total(numFile) = length(validCells);

        infoEuler3D{numFile, 1} = array2table(vertcat(meanNumNeighPerSurfaceRealization(:, 1:minNumberOfSurfaceRatios), stdNumNeighPerSurfaceRealization(:, 1:minNumberOfSurfaceRatios), mean_PercScutoids(:, 1:minNumberOfSurfaceRatios), std_PercScutoids(:, 1:minNumberOfSurfaceRatios), mean_apicoBasalTransitions(:, 1:minNumberOfSurfaceRatios), std_apicoBasalTransitions(:, 1:minNumberOfSurfaceRatios), surfaceRatioOfGland(:, 1:minNumberOfSurfaceRatios))','VariableNames',{'mean_neigh3D','std_neigh3D','mean_PercScutoids','std_PercScutoids', 'mean_apicoBasalTransitions', 'std_apicoBasalTransitions','Surface_Ratio'});
        numNeighPerSurface{numFile, 1} = array2table(numNeighPerSurfaceRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
        numNeighAccumPerSurfaces{numFile, 1} = array2table(numNeighAccumPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
        numNeighOfNeighPerSurface{numFile, 1} = array2table(numNeighOfNeighPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
        numNeighOfNeighAccumPerSurface{numFile, 1} = array2table(numNeighOfNeighAccumPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
        areaCellsPerSurface{numFile, 1} = array2table(areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios) ./ mean(areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios)),'VariableNames',namesSR);
        volumePerSurface{numFile, 1} = array2table(volumePerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios) ./ mean(volumePerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios)),'VariableNames',namesSR);

        %columns = {'NumCell','numberOfNeighbours','area','numberOfNeighbours','area','numberOfNeighbours','area','numberOfNeighbours','area','numberOfNeighbours','area','numberOfNeighbours','area','numberOfNeighbours','area','numberOfNeighbours','area','numberOfNeighbours','area','numberOfNeighbours','area','numberOfNeighbours','area'};
        %surfaceRatioOfGland_Text = {'Surface ratio' surfaceRatioOfGland(1) [] surfaceRatioOfGland(2) [] surfaceRatioOfGland(3)  [] surfaceRatioOfGland(4)  [] surfaceRatioOfGland(5)  [] surfaceRatioOfGland(6)  [] surfaceRatioOfGland(7)  [] surfaceRatioOfGland(8)  [] surfaceRatioOfGland(9)  [] surfaceRatioOfGland(10)  [] surfaceRatioOfGland(11) []};
        totalDataBasal = [totalDataBasal; validCells', areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios) ./ mean(areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios)), numNeighPerSurfaceRealization(validCells, 1:minNumberOfSurfaceRatios)];
        %T = totalDataBasal(:, [1 13 2 14 3 15 4 16 5 17 6 18 7 19 8 20 9 21 10 22 11 23 12]);
        %xlswrite([num2str(numFile) '_tableAreaSidesBasal.xls'], [surfaceRatioOfGland_Text; columns; num2cell(T)])
        totalDataAccum = [totalDataAccum; validCells', areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios) ./ mean(areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios)), numNeighAccumPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios)];
        %T = totalDataAccum(:, [1 13 2 14 3 15 4 16 5 17 6 18 7 19 8 20 9 21 10 22 11 23 12]);
        %xlswrite([num2str(numFile) '_tableAreaSidesAccum.xls'], [surfaceRatioOfGland_Text; columns; num2cell(T)])

        realSRBasal = [realSRBasal; validCells', areaCellsPerSurfaceRealization(validCells,[1 end]) ./ mean(areaCellsPerSurfaceRealization(validCells, [1 end])), numNeighPerSurfaceRealization(validCells, [1 end])];
        realSRAccum = [realSRAccum; validCells', areaCellsPerSurfaceRealization(validCells,[1 end]) ./ mean(areaCellsPerSurfaceRealization(validCells,[1 end])), numNeighAccumPerSurfacesRealization(validCells, [1 end])];
        %Scutoids per number of sides
        [meanWinningPerSidePerFile{numFile, 1}, cellsPerSide{numFile}] = calculateMeanWinning3DNeighbours(numNeighAccumPerSurfacesRealization(:, 1:minNumberOfSurfaceRatios), validCells, minNumberOfSurfaceRatios);
        clearvars 'meanNeighsScutoidsPerSF_ValidCells' 'neighbours'
    end
    
    %%
    meanNumNeighPerSurfaceRealization = mean(numNeighAccumPerSurfacesRealization(validCells, :), 1);
    numCells = repmat(length(validCells), 1, size(meanNumNeighPerSurfaceRealization, 2));
    stdNumNeighPerSurfaceRealization = std(numNeighAccumPerSurfacesRealization(validCells, :), 1);
    totalAreaPerSR = sum(areaCellsPerSurfaceRealization(validCells, :));
    
    mean_PercScutoids = mean(percentageScutoids(validCells, :), 1);
    std_PercScutoids = std(percentageScutoids(validCells, :), 1);
    
    mean_apicoBasalTransitions = mean(apicoBasalTransitions(validCells, :), 1);
    std_apicoBasalTransitions = std(apicoBasalTransitions(validCells, :), 1);
    
    surfaceRatioOfGland_real = vertcat(infoPerSurfaceRatio{:, 7})'; 
    totalPartitions = 10;
    initialPartitions = (1:(totalPartitions-1))/totalPartitions;
    surfaceRatioOfGland = surfaceRatioOfGland_real;
    surfaceRatioOfGland(2:10) = initialPartitions * (surfaceRatioOfGland_real(end) - 1) + 1;
    
    
    infoEuler3D{numFile, 1} = array2table(vertcat(meanNumNeighPerSurfaceRealization, stdNumNeighPerSurfaceRealization, mean_PercScutoids, std_PercScutoids, mean_apicoBasalTransitions, std_apicoBasalTransitions, surfaceRatioOfGland)','VariableNames',{'mean_neigh3D','std_neigh3D','mean_PercScutoids','std_PercScutoids', 'mean_apicoBasalTransitions', 'std_apicoBasalTransitions','Surface_Ratio'});
    numNeighPerSurface{numFile, 1} = array2table(numNeighPerSurfaceRealization(validCells, :),'VariableNames',namesSR);
    numNeighAccumPerSurfaces{numFile, 1} = array2table(numNeighAccumPerSurfacesRealization(validCells, :),'VariableNames',namesSR);
    numNeighOfNeighPerSurface{numFile, 1} = array2table(numNeighOfNeighPerSurfacesRealization(validCells, :),'VariableNames',namesSR);
    numNeighOfNeighAccumPerSurface{numFile, 1} = array2table(numNeighOfNeighAccumPerSurfacesRealization(validCells, :),'VariableNames',namesSR);
    areaCellsPerSurface{numFile, 1} = array2table(areaCellsPerSurfaceRealization(validCells,:),'VariableNames',namesSR);
    volumePerSurface{numFile, 1} = array2table(volumePerSurfaceRealization(validCells,:),'VariableNames',namesSR);
    
    %Scutoids per number of sides
    [meanWinningPerSidePerFile{numFile, 1}, cellsPerSide{numFile}] = calculateMeanWinning3DNeighbours(numNeighAccumPerSurfacesRealization, validCells);
    clearvars 'meanNeighsScutoidsPerSF_ValidCells' 'neighbours'

%     save('Results/salivaryGland_Info_20_05_2019.mat', 'mean_PercScutoids_basal', 'infoEuler3D', 'numNeighPerSurface', 'numNeighAccumPerSurfaces', 'numNeighOfNeighPerSurface', 'numNeighOfNeighAccumPerSurface', 'areaCellsPerSurface', 'volumePerSurface', 'numCells_Total', 'neighbours_apical', 'neighbours_basal', 'neighbours_total', 'polygon_distribution_apical', 'polygon_distribution_basal','meanVolumeMicronsPerGland', 'mean_apicoBasalTransitionsPerGland', 'mean_neighsAccum');
% else
%     load('Results/salivaryGland_Info_20_05_2019.mat');
% end

nTotal_apicobasalT = [mean_apicoBasalTransitionsPerGland{:}; mean_neighsAccum{:}; mean_Scutoids{:}; surfaceRatiosPerFile{:}]';
%% NTotal vs ApicoBasal
figure; 
%hold on;
for numPoint = 1:size(nTotal_apicobasalT, 1)
    plot(nTotal_apicobasalT(numPoint, 1), nTotal_apicobasalT(numPoint, 2), 'o', 'MarkerSize', 5, 'Color', [151 238 152]/255, 'MarkerFaceColor', [151 238 152]/255, 'MarkerEdgeColor', [151 238 152]/255)
    hold on;
end
% %% NTotal vs Scutoids
% figure; 
% for numPoint = 1:size(nTotal_apicobasalT, 1)
%     plot(nTotal_apicobasalT(numPoint, 3), nTotal_apicobasalT(numPoint, 2), 'ko','MarkerSize',5, 'MarkerFaceColor',[151 238 152]/255)
%     hold on;
% end
% xlabel('% Scutoids')
% ylabel('N Total')
% legend('S. Gland');
% %% NTotal vs SR
% figure; 
% for numPoint = 1:size(nTotal_apicobasalT, 1)
%     plot(nTotal_apicobasalT(numPoint, 4), nTotal_apicobasalT(numPoint, 2), 'ko','MarkerSize',5, 'MarkerFaceColor',[151 238 152]/255, )
%     hold on;
% end
% xlabel('Surface ratio')
% ylabel('N Total')
% legend('S. Gland');
%%
mean_apicoBasalTransitions_final
dim = ndims(meanWinningPerSidePerFile{1});          %# Get the number of dimensions for your arrays
M = cat(dim+1,meanWinningPerSidePerFile{:});        %# Convert to a (dim+1)-dimensional matrix
meanWinningPerSide_Total = mean(M,dim+1, 'omitnan');  %# Get the mean across arrays

infoEuler3DCat = cat(1, infoEuler3D{:,1});

%figure;
myfittypeLn=fittype('6 + b*log(x)',...
'dependent', {'y'}, 'independent',{'x'},...
'coefficients', {'b'});
myfittypePoly=fittype('6 +b*x',...
'dependent', {'y'}, 'independent',{'x'},...
'coefficients', {'b'});
goodnessln = cell(size(infoEuler3D, 1), 1);
predD = cell(size(infoEuler3D, 1), 1);
outputlog = cell(size(infoEuler3D, 1), 1);
rSquaresln = zeros(size(infoEuler3D, 1),1);
% coefAlog = zeros(size(infoEuler3D, 1),1);
coefBlog = zeros(size(infoEuler3D, 1),1);
goodnessPol = cell(size(infoEuler3D, 1), 1);
outputPol = cell(size(infoEuler3D, 1), 1);
rSquaresPol = zeros(size(infoEuler3D, 1),1);
% coefAPol = zeros(size(infoEuler3D, 1),1);
coefBPol = zeros(size(infoEuler3D, 1),1);
goodnessPol_fromLn = cell(size(infoEuler3D, 1), 1);
outputPol_fromLn = cell(size(infoEuler3D, 1), 1);
rSquaresPol_fromLn = zeros(size(infoEuler3D, 1),1);

for numPoint = 1:size(infoEuler3D, 1)
    
    infoEulerActual = infoEuler3D{numPoint,1};
    if isempty(infoEulerActual)
        continue;
    end
    %infoEulerActual = infoEulerActual([1 6 11], :);
    figure;
    errorbar(infoEulerActual.Surface_Ratio,infoEulerActual.mean_neigh3D,infoEulerActual.std_neigh3D,'-o','MarkerSize',5,...
        'MarkerEdgeColor','black','MarkerFaceColor','blue');
    
    [myfitLn, goodnessln{numPoint}, outputlog{numPoint}] = fit(infoEulerActual.Surface_Ratio,infoEulerActual.mean_neigh3D,myfittypeLn,'StartPoint',6);
    predD{numPoint} = predint(myfitLn,infoEulerActual.Surface_Ratio,0.95,'functional','on');
    
    rSquaresln(numPoint) = goodnessln{numPoint}.rsquare;
%     coefAlog(numPoint) = myfitLn.a;
    coefBlog(numPoint) = myfitLn.b;
    
    [myfitPol_fromLog, goodnessPol_fromLn{numPoint}, outputPol_fromLn{numPoint}] = fit(infoEulerActual.Surface_Ratio, myfitLn(infoEulerActual.Surface_Ratio), myfittypePoly,'StartPoint',6); 

    rSquaresPol_fromLn(numPoint) = goodnessPol_fromLn{numPoint}.rsquare;
    
    hold on; plot(myfitLn);
    title('euler neighbours 3D')
    xlabel('surface ratio')
    ylabel('neighbours total')
    xlim([1, 8]);
    ylim([0,15]);
    hold off;
    [myfitPol, goodnessPol{numPoint}, outputPol{numPoint}] = fit(infoEulerActual.Surface_Ratio, infoEulerActual.mean_neigh3D, myfittypePoly,'StartPoint',6);
    
    rSquaresPol(numPoint) = goodnessPol{numPoint}.rsquare;
%     coefAPol(numPoint) = myfitPol.a;
    coefBPol(numPoint) = myfitPol.b;
    
    figure(100)
    hold on
    plot(infoEulerActual.mean_PercScutoids,infoEulerActual.mean_neigh3D,'o','MarkerSize',5,...
        'MarkerEdgeColor','black','MarkerFaceColor','blue');
    xlabel('% Scutoids')
    ylabel('TotalNeighs')
    hold off
    
    figure(101)
    hold on
    plot(infoEulerActual.mean_apicoBasalTransitions, infoEulerActual.mean_neigh3D,'o','MarkerSize',5,...
        'MarkerEdgeColor','black','MarkerFaceColor',[151 238 152]/255);
    xlabel('Apico-Basal transitions')
    ylabel('TotalNeighs')
    hold off
    
    figure(102)
    hold on
    plot(infoEulerActual.mean_PercScutoids,infoEulerActual.mean_apicoBasalTransitions,'o','MarkerSize',5,...
        'MarkerEdgeColor','black','MarkerFaceColor','blue');
    xlabel('% Scutoids')
    ylabel('Apico-Basal transitions')
    hold off
end

% meanRsquareLn = mean(rSquaresln)
% stdRsquareLn = std(rSquaresln);
% meanCoefALn = mean(coefAlog);
% stdCoefALn = std(coefAlog);
% meanCoefBLn = mean(coefBlog);
% stdCoefBLn = std(coefBlog);
% 
% meanRsquarePol = mean(rSquaresPol)
% stdRsquarePol = std(rSquaresPol);
% meanCoefAPol = mean(coefAPol);
% stdCoefAPol = std(coefAPol);
% meanCoefBPol = mean(coefBPol);
% stdCoefBPol = std(coefBPol);
% 
% meanRsquarePol_fromLn = mean(rSquaresPol_fromLn)

% myfittypeLn=fittype('6 +b*log(x)',...
% 'dependent', {'y'}, 'independent',{'x'},...
% 'coefficients', {'b'});

% a = mean(coefAlog);
% b = mean(coefBlog);
% actualFit = cfit(myfittypeLn, a, b);

% % a = mean(coefAPol);
% % b = mean(coefBPol);
% % actualFit = cfit(myfittypePoly, a, b);
% 
% % figure; 
% allPoints = vertcat(infoEuler3D{:, 1});
% 
% for numPoint = 1:size(allPoints, 1)
%     hold on;
%     plot(allPoints.Surface_Ratio(numPoint) , allPoints.mean_neigh3D(numPoint), 'o', 'Color', [151 238 152]/255, 'MarkerFaceColor', [151 238 152]/255, 'LineWidth', 0.2, 'MarkerSize', 5);
% end
% %plot(actualFit, [1 11], [6 7]); %LineWidth 2
% plot(actualFit, [0 1], [0.75 1]);
% ylim([0 15])
% xlim([0 15])
% x = [0 16];
% y = [6 6];
% line(x, y, 'Color', 'red', 'LineStyle', '--')
% title(num2str(mean(meanRsquareLn)))
% xlabel('surface ratio')
% ylabel('neighbours total')
% figure;
% 
% %x = ((0:19) * 11) + [1 6 11]';
% for numPoint = 1:size(infoEuler3DCat, 1)%x(:)
%     
%     hold on;
%     plot(infoEuler3DCat.Surface_Ratio(numPoint), infoEuler3DCat.mean_neigh3D(numPoint), 'o', 'Color', [151 238 152]/255, 'MarkerFaceColor', [151 238 152]/255, 'LineWidth', 0.2, 'MarkerSize', 5)
% end
% title('euler neighbours 3D')
% xlabel('surface ratio')
% ylabel('neighbours total')
% xlim([1, 15]);
% ylim([0,15]);
% 
% 
% 
% [myfitLn, goodnessln, outputln]=fit(infoEuler3DCat.Surface_Ratio, infoEuler3DCat.mean_neigh3D, myfittypeLn,'StartPoint',6);
% 
% yAxis = [unique(infoEuler3DCat.Surface_Ratio)' max(unique(infoEuler3DCat.Surface_Ratio)):11 11];
% preD = predint(myfitLog10, yAxis,0.95, 'observation','off');
% plot(yAxis,preD,'--','Color',[0.00,0.80,0.00])
% hold on; plot(myfitLog10, [1 11], [6 myfitLog10(11)]);
% 
% 
% [myfitPol, goodnessPol, outputPol]=fit(infoEuler3DCat.Surface_Ratio, infoEuler3DCat.mean_neigh3D,myfittypePoly,'StartPoint',[6, 1]);
% hold on; plot(myfitPol);

getStatsAndRepresentationsEulerLewis3D(numNeighOfNeighPerSurface,numNeighOfNeighAccumPerSurface,numNeighPerSurface,numNeighAccumPerSurfaces,areaCellsPerSurface,volumePerSurface,'Results\SalivaryGlands\', surfaceRatiosExtrapolatedFrom3D, 0);