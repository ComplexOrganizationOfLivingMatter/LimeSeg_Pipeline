
% files = dir('**/Salivary gland/**/Results/dividedGland/glandDividedInSurfaceRatios.mat');
% for numFile = 1:length(files)
%     load(fullfile(files(numFile).folder, 'glandDividedInSurfaceRatios.mat'))
%     %figure; paint3D( ismember(imageOfSurfaceRatios{numPartition, 3}, validCells) .* imageOfSurfaceRatios{numPartition, 3}, [], colours);
%     intermidScutoids = sum(infoPerSurfaceRatio{6, 4}.Scutoids) / size(infoPerSurfaceRatio{6, 4}, 1);
%     finalScutoids = sum(infoPerSurfaceRatio{11, 4}.Scutoids) / size(infoPerSurfaceRatio{11, 4}, 1);
%     logScutoids(numFile) = intermidScutoids/finalScutoids;
% end


clear all
close all

addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files = dir('**/Salivary gland/**/Results/3d_layers_info.mat');

nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
files = files(nonDiscardedFiles);

%resultsFileName = '3d_layers_info.mat';
resultsFileName = 'glandDividedInSurfaceRatios.mat';
%resultsFileName = 'glandDividedInSurfaceRatios_PredefinedSR.mat';

for numFile = 1:length(files)
    files(numFile).folder
    load(fullfile(files(numFile).folder, 'valid_cells'), 'validCells', 'noValidCells');
    load(fullfile(files(numFile).folder, 'apical', 'verticesInfo'), 'cylindre2DImage');
    load(fullfile(files(numFile).folder, '3d_layers_info'), 'colours');
%     areaApical = sum(sum(ismember(cylindre2DImage, validCells)));
    
    load(fullfile(files(numFile).folder, 'basal', 'verticesInfo'), 'cylindre2DImage');
%     areaBasal = sum(sum(ismember(cylindre2DImage, validCells)));  
%     surfaceRatioGlans(numFile) = areaBasal / areaApical;
    
    if exist(fullfile(files(numFile).folder, 'dividedGland' ,'glandDividedInSurfaceRatios.mat'), 'file') > 0
        load(fullfile(files(numFile).folder, 'dividedGland', 'glandDividedInSurfaceRatios.mat'))
        imageOfSurfaceRatios = infoPerSurfaceRatio;
        parfor numPartition = 2:10
            unrollTube(imageOfSurfaceRatios{numPartition, 3}, fullfile(files(numFile).folder, 'dividedGland', ['gland_SR_' num2str(imageOfSurfaceRatios{numPartition, 2})]), noValidCells, colours, 1);
        end
    end
end
% mean(surfaceRatioGlans)