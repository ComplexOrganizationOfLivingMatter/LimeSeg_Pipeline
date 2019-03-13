addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))
%addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files = dir('**/Salivary gland/**/Results/3d_layers_info.mat');

% for numFile = 1:length(files)
%     if contains(lower(files(numFile).folder), 'discarded') == 0
%         selpath = files(numFile).folder;
%         divideObjectInSurfaceRatios(selpath);
%     end
% end

resultsFileName = 'glandDividedInSurfaceRatios.mat';
%resultsFileName = 'glandDividedInSurfaceRatios_PredefinedSR.mat';

for numFile = 1:length(files)
    if exist(fullfile(files(numFile).folder, 'dividedGland' ,'glandDividedInSurfaceRatios.mat'), 'file') > 0
        
        files(numFile).folder
        
        load(fullfile(files(numFile).folder, 'valid_cells'));
        filesOf2DUnroll = dir(fullfile(files(numFile).folder, '**', 'verticesInfo.mat'));
        if length(filesOf2DUnroll) == 11
            load(fullfile(filesOf2DUnroll(1).folder, 'verticesInfo.mat'));
            apicalLayer = cylindre2DImage;
            neighboursApical = arrayfun(@(x) unique(newVerticesNeighs2D(any(ismember(newVerticesNeighs2D, x), 2), :)), 1:max(newVerticesNeighs2D(:)), 'UniformOutput', false);
            neighboursApical = cellfun(@(x, y) x(x ~= y), neighboursApical, num2cell(1:length(neighboursApical)), 'UniformOutput', false);
            neighboursApical = cellfun(@(x) x(:), neighboursApical, 'UniformOutput', false);
            
            load(fullfile(filesOf2DUnroll(2).folder, 'verticesInfo.mat'));
            basalLayer = cylindre2DImage;
            neighboursBasal = arrayfun(@(x) unique(newVerticesNeighs2D(any(ismember(newVerticesNeighs2D, x), 2), :)), 1:max(newVerticesNeighs2D(:)), 'UniformOutput', false);
            neighboursBasal = cellfun(@(x, y) x(x ~= y), neighboursBasal, num2cell(1:length(neighboursBasal)), 'UniformOutput', false);
            neighboursBasal = cellfun(@(x) x(:), neighboursBasal, 'UniformOutput', false);
            
            for numSR = 3:11
                numSR
                load(fullfile(filesOf2DUnroll(numSR).folder, 'final3DImg.mat'), 'img3d');
                load(fullfile(filesOf2DUnroll(numSR).folder, 'verticesInfo.mat'), 'cylindre2DImage', 'newVerticesNeighs2D');
                midLayer = cylindre2DImage;
                neighboursMid = arrayfun(@(x) unique(newVerticesNeighs2D(any(ismember(newVerticesNeighs2D, x), 2), :)), 1:max(newVerticesNeighs2D(:)), 'UniformOutput', false);
                neighboursMid = cellfun(@(x, y) x(x ~= y), neighboursMid, num2cell(1:length(neighboursMid)), 'UniformOutput', false);
                neighboursMid = cellfun(@(x) x(:), neighboursMid, 'UniformOutput', false);
            
                neighbours_data = table(neighboursApical, neighboursMid);
                neighbours_data.Properties.VariableNames = {'Apical','Basal'};
                [cellularFeatures_ApicalToBasal{numSR - 2}, meanSurfaceRatio(numSR-2)] = calculate_CellularFeatures(neighbours_data, neighboursApical, neighboursMid, apicalLayer, midLayer, img3d, noValidCells, validCells, [], []);
                
%                 neighbours_data = table(neighboursMid, neighboursApical);
%                 neighbours_data.Properties.VariableNames = {'Apical','Basal'};
%                 [cellularFeatures_BasalToApical{numSR - 2}, ~] = calculate_CellularFeatures(neighbours_data, neighboursMid, neighboursBasal, midLayer, basalLayer, img3d, noValidCells, validCells, [], []);
            end
            
            save(fullfile(files(numFile).folder, 'glandDividedInSurfaceRatios_AllUnrollFeatures.mat'), 'cellularFeatures_BasalToApical', 'cellularFeatures_ApicalToBasal', 'meanSurfaceRatio');
            
            save(fullfile(files(numFile).folder, 'glandDividedInSurfaceRatios_AllUnrollFeatures.mat'), 'cellularFeatures_ApicalToBasal', 'meanSurfaceRatio');
        else
%             load(fullfile(files(numFile).folder, '3d_layers_info'), 'colours');
%             load(fullfile(files(numFile).folder, 'dividedGland', 'glandDividedInSurfaceRatios.mat'));
%             
%             imageOfSurfaceRatios = infoPerSurfaceRatio;
%             parfor numPartition = 2:10
%                 unrollTube(imageOfSurfaceRatios{numPartition, 3}, fullfile(files(numFile).folder, 'dividedGland', ['gland_SR_' num2str(imageOfSurfaceRatios{numPartition, 2})]), noValidCells, colours, 1);
%             end
        end
    end
end