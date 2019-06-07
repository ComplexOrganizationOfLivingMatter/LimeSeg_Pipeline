
clear all
close all

addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files = dir('**/Salivary gland/**/Results/3d_layers_info.mat');

nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0 && contains(lower(x), 'wildtype'), {files.folder});
files = files(nonDiscardedFiles);

%resultsFileName = '3d_layers_info.mat';
%resultsFileName = 'glandDividedInSurfaceRatios.mat';
%resultsFileName = 'glandDividedInSurfaceRatios_PredefinedSR.mat';
resultsFileName = 'glandDividedInSurfaceRatios_AllUnrollFeatures.mat';

numberOfSurfaceRatios = 11;
namesSR = arrayfun(@(x) ['sr' strrep(num2str(x),'.','_')],1:numberOfSurfaceRatios,'UniformOutput', false);
%namesSR = {'sr1' 'sr2'};

steps = 2.5/(minNumberOfSurfaceRatios-1);
surfaceRatiosExtrapolatedFrom3D = 1:steps:((steps*(minNumberOfSurfaceRatios-1))+1);

totalDataBasal = [];
totalDataAccum = [];

realSRBasal = [];
realSRAccum = [];
meanVolumeMicronsPerGland = zeros(length(files),2);


if ~exist('Results/salivaryGland_Info_20_05_2019.mat','file')
    for numFile = 1:length(files)
        files(numFile).folder

        load(fullfile(files(numFile).folder, files(numFile).name));
        load(fullfile(files(numFile).folder, 'valid_cells.mat'));
        if exist(fullfile(files(numFile).folder, resultsFileName), 'file')
            load(fullfile(files(numFile).folder, resultsFileName))
        else
            continue
            %divideObjectInSurfaceRatios(files(numFile).folder);
        end
        load([files(numFile).folder '\unrolledGlands\gland_SR_basal\final3DImg.mat'],'img3dComplete')
        
        volume3d = regionprops3(img3dComplete,'Volume');
        volume3d = cat(1,volume3d.Volume);
        meanVolumeMicronsPerGland(numFile,1) = mean(volume3d(validCells)); 
        meanVolumeMicronsPerGland(numFile,2) = std(volume3d(validCells)); 
        if exist('infoPerSurfaceRatio', 'var')
            numberOfSurfaceRatios = size(infoPerSurfaceRatio, 1);
        else
            numberOfSurfaceRatios = 2;
        end

        meanSR(numFile, 1) = infoPerSurfaceRatio{end, 2};
        meanSR(numFile, 2) = infoPerSurfaceRatio{end, 7};
        meanSR(numFile, 3) = infoPerSurfaceRatio{minNumberOfSurfaceRatios, 2};
        meanSR(numFile, 4) = infoPerSurfaceRatio{minNumberOfSurfaceRatios, 7};
        if exist('neighboursOfAllSurfaces', 'var') == 0
            neighboursOfAllSurfaces = cell(numberOfSurfaceRatios, 1);
            filesOf2DUnroll = dir(fullfile(files(numFile).folder, '**', 'verticesInfo.mat'));
            for numSR = 1:numberOfSurfaceRatios
                load(fullfile(filesOf2DUnroll(numSR).folder, 'verticesInfo.mat'), 'newVerticesNeighs2D');

                if numSR == 1
                    idToSave = 1;
                elseif numSR == 2
                    idToSave = size(infoPerSurfaceRatio, 1);
                else
                    idToSave = numSR - 1;
                end

                neighboursOfAllSurfaces{idToSave} = getNeighboursFromVertices(newVerticesNeighs2D);
            end
            
            neighboursOfAllSurfaces{idToSave} = getNeighboursFromVertices(newVerticesNeighs2D);
        end
    end
    
    neighsSurface = cell(numberOfSurfaceRatios,1);
    neighsAccumSurfaces = cell(numberOfSurfaceRatios,1);
    percentageScutoids = cell(numberOfSurfaceRatios, 1);
    apicoBasalTransitions = cell(numberOfSurfaceRatios, 1);
    areaCells = cell(numberOfSurfaceRatios,1);
    volumes = cell(numberOfSurfaceRatios,1);
    
    neighsSurface{1} = neighboursOfAllSurfaces{1};
    neighsAccumSurfaces{1} = neighboursOfAllSurfaces{1};
    percentageScutoids{1} = cellfun(@(x, y) ~isequal(x,y), neighsSurface{1}, neighsAccumSurfaces{1});
    apicoBasalTransitions{1} = cellfun(@(x, y) length(setxor(x,y)), neighsSurface{1}, neighsSurface{1});
    
    infoOfCells = infoPerSurfaceRatio{1, 4};
    infoOfCells = infoOfCells{:};
    
    areaCells(1) = {infoOfCells.Basal_area};
    volumes(1) = {infoOfCells.Volume};
    
    for idSR = 2:numberOfSurfaceRatios
        neighsSurface{idSR} = neighboursOfAllSurfaces{idSR};
        neighsAccumSurfaces{idSR} = cellfun(@(x,y) unique([x;y]),neighsAccumSurfaces{idSR-1},neighsSurface{idSR},'UniformOutput',false);
        percentageScutoids{idSR} = cellfun(@(x, y) ~isempty(setxor(x,y)), neighsSurface{1}, neighsSurface{idSR});
        apicoBasalTransitions{idSR} = cellfun(@(x, y) length(setxor(x,y)), neighsSurface{1}, neighsSurface{idSR});
        
        infoOfCells = infoPerSurfaceRatio{idSR, 4};
        infoOfCells = infoOfCells{:};
        areaCells(idSR) = {infoOfCells.Basal_area};
        volumes(idSR) = {infoOfCells.Volume};
    end
    
    areaCellsPerSurfaceRealization = cat(2,areaCells{:});
    volumePerSurfaceRealization = cat(2,volumes{:});
    neighsSurface = cat(1,neighsSurface{:})';
    neighsAccumSurfaces = cat(1,neighsAccumSurfaces{:})';
    percentageScutoids = cat(1,percentageScutoids{:})';
    apicoBasalTransitions = cat(1,apicoBasalTransitions{:})';
    
    numNeighPerSurfaceRealization = cellfun(@(x) length(x),neighsSurface);
    numNeighAccumPerSurfacesRealization = cellfun(@(x) length(x),neighsAccumSurfaces);
    
    numNeighOfNeighPerSurfacesRealization = zeros(size(neighsSurface));
    numNeighOfNeighAccumPerSurfacesRealization = zeros(size(neighsSurface));
    for nSR = 1:numberOfSurfaceRatios
        numNeighOfNeighPerSurfacesRealization(:,nSR) = cellfun(@(x) sum(vertcat(numNeighPerSurfaceRealization(x,nSR)))/length(x),neighsSurface(:,nSR));
        numNeighOfNeighAccumPerSurfacesRealization(:,nSR) = cellfun(@(x) sum(vertcat(numNeighAccumPerSurfacesRealization(x,nSR)))/length(x),neighsAccumSurfaces(:,nSR));
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

    save('Results/salivaryGland_Info_20_05_2019.mat', 'mean_PercScutoids_basal', 'infoEuler3D', 'numNeighPerSurface', 'numNeighAccumPerSurfaces', 'numNeighOfNeighPerSurface', 'numNeighOfNeighAccumPerSurface', 'areaCellsPerSurface', 'volumePerSurface', 'numCells_Total', 'neighbours_apical', 'neighbours_basal', 'neighbours_total', 'polygon_distribution_apical', 'polygon_distribution_basal','meanVolumeMicronsPerGland', 'mean_apicoBasalTransitionsPerGland', 'mean_neighsAccum');
else
    load('Results/salivaryGland_Info_20_05_2019.mat');
end

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
myfittypeLog10=fittype('a +b*log10(x)',...
'dependent', {'y'}, 'independent',{'x'},...
'coefficients', {'a','b'});
myfittypePoly=fittype('a +b*x',...
'dependent', {'y'}, 'independent',{'x'},...
'coefficients', {'a','b'});
goodnesslog = cell(size(infoEuler3D, 1), 1);
outputlog = cell(size(infoEuler3D, 1), 1);
rSquareslog = zeros(size(infoEuler3D, 1),1);
coefAlog = zeros(size(infoEuler3D, 1),1);
coefBlog = zeros(size(infoEuler3D, 1),1);
goodnessPol = cell(size(infoEuler3D, 1), 1);
outputPol = cell(size(infoEuler3D, 1), 1);
rSquaresPol = zeros(size(infoEuler3D, 1),1);
coefAPol = zeros(size(infoEuler3D, 1),1);
coefBPol = zeros(size(infoEuler3D, 1),1);
goodnessPol_fromLog = cell(size(infoEuler3D, 1), 1);
outputPol_fromLog = cell(size(infoEuler3D, 1), 1);
rSquaresPol_fromLog = zeros(size(infoEuler3D, 1),1);

for numPoint = 1:size(infoEuler3D, 1)
    
    infoEulerActual = infoEuler3D{numPoint,1};
    if isempty(infoEulerActual)
        continue;
    end
    %infoEulerActual = infoEulerActual([1 6 11], :);
    figure;
    errorbar(infoEulerActual.Surface_Ratio,infoEulerActual.mean_neigh3D,infoEulerActual.std_neigh3D,'-o','MarkerSize',5,...
        'MarkerEdgeColor','black','MarkerFaceColor','blue');
    
    [myfitLog10, goodnesslog{numPoint}, outputlog{numPoint}] = fit(infoEulerActual.Surface_Ratio,infoEulerActual.mean_neigh3D,myfittypeLog10,'StartPoint',[6,1]);
    
    rSquareslog(numPoint) = goodnesslog{numPoint}.rsquare;
    coefAlog(numPoint) = myfitLog10.a;
    coefBlog(numPoint) = myfitLog10.b;
    
    [myfitPol_fromLog, goodnessPol_fromLog{numPoint}, outputPol_fromLog{numPoint}] = fit(infoEulerActual.Surface_Ratio, myfitLog10(infoEulerActual.Surface_Ratio), myfittypePoly,'StartPoint',[6,1]); 
    
    rSquaresPol_fromLog(numPoint) = goodnessPol_fromLog{numPoint}.rsquare;
    
    hold on; plot(myfitLog10);
    title('euler neighbours 3D')
    xlabel('surface ratio')
    ylabel('neighbours total')
    xlim([1, 8]);
    ylim([0,15]);
    hold off;
    [myfitPol, goodnessPol{numPoint}, outputPol{numPoint}] = fit(infoEulerActual.Surface_Ratio, infoEulerActual.mean_neigh3D, myfittypePoly,'StartPoint',[6,1]);
    
    rSquaresPol(numPoint) = goodnessPol{numPoint}.rsquare;
    coefAPol(numPoint) = myfitPol.a;
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

meanRsquareLog = mean(rSquareslog)
stdRsquareLog = std(rSquareslog);
meanCoefALog = mean(coefAlog);
stdCoefALog = std(coefAlog);
meanCoefBLog = mean(coefBlog);
stdCoefBLog = std(coefBlog);

meanRsquarePol = mean(rSquaresPol)
stdRsquarePol = std(rSquaresPol);
meanCoefAPol = mean(coefAPol);
stdCoefAPol = std(coefAPol);
meanCoefBPol = mean(coefBPol);
stdCoefBPol = std(coefBPol);

meanRsquarePol_fromLog = mean(rSquaresPol_fromLog)

myfittypeLog10=fittype('a +b*log10(x)',...
'dependent', {'y'}, 'independent',{'x'},...
'coefficients', {'a', 'b'});

a = mean(coefAlog);
b = mean(coefBlog);
actualFit = cfit(myfittypeLog10, a, b);


a = mean(coefAPol);
b = mean(coefBPol);
actualFit = cfit(myfittypePoly, a, b);

figure; 
allPoints = vertcat(infoEuler3D{:, 1});
for numPoint = 1:size(allPoints, 1)
    hold on;
    plot(allPoints.Surface_Ratio(numPoint) , allPoints.mean_neigh3D(numPoint), 'o', 'Color', [151 238 152]/255, 'LineWidth', 2, 'MarkerSize', 5);
end
%plot(actualFit, [1 11], [6 7]); %LineWidth 2
plot(actualFit, [0 1], [0.75 1]);
ylim([0 15])
xlim([0 15])
x = [0 16];
y = [6 6];
line(x, y, 'Color', 'red', 'LineStyle', '--')
title(num2str(mean(meanRsquareLog)))
xlabel('surface ratio')
ylabel('neighbours total')
figure;


for numPoint = 1:size(infoEuler3DCat, 1)
    hold on;
    plot(infoEuler3DCat(numPoint, 3), infoEuler3DCat(numPoint, 1), '*k')
end
title('euler neighbours 3D')
xlabel('surface ratio')
ylabel('neighbours total')
xlim([1, 15]);
ylim([0,15]);



[myfitLog10, goodnesslog, outputlog]=fit(infoEuler3DCat.Surface_Ratio, infoEuler3DCat.mean_neigh3D, myfittypeLog10,'StartPoint',[6, 1]);
hold on; plot(myfitLog10);
[myfitPol, goodnessPol, outputPol]=fit(infoEuler3DCat.Surface_Ratio, infoEuler3DCat.mean_neigh3D,myfittypePoly,'StartPoint',[6, 1]);
hold on; plot(myfitPol);


getStatsAndRepresentationsEulerLewis3D(numNeighOfNeighPerSurface,numNeighOfNeighAccumPerSurface,numNeighPerSurface,numNeighAccumPerSurfaces,areaCellsPerSurface,volumePerSurface,'Results/SalivaryGlands/', 1:numberOfSurfaceRatios);