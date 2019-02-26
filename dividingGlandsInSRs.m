addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))
%addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files = dir('**/Salivary gland/**/Results/3d_layers_info.mat');

for numFile = 1:length(files)
    if contains(lower(files(numFile).folder), 'discarded') == 0
        selpath = files(numFile).folder;
        load(fullfile(selpath, '3d_layers_info.mat'));
        load(fullfile(selpath, 'valid_cells.mat'));
        
        load(fullfile(selpath, 'apical', 'verticesInfo.mat'), 'newVerticesNeighs2D');
        apicalNeighs = newVerticesNeighs2D;
        apicalRealNeighs = arrayfun(@(x) sum(any(ismember(apicalNeighs, x), 2)), 1:max(validCells));
        
        load(fullfile(selpath, 'basal', 'verticesInfo.mat'), 'newVerticesNeighs2D');
        basalNeighs = newVerticesNeighs2D;
        basalRealNeighs = arrayfun(@(x) sum(any(ismember(basalNeighs, x), 2)), 1:max(validCells));
        
        [infoPerSurfaceRatio, neighbours] = divideObjectInSurfaceRatios(labelledImage, basalLayer, apicalLayer, validCells, noValidCells, colours, selpath, basalRealNeighs, apicalRealNeighs);

        save(fullfile(selpath, 'dividedGland', 'glandDividedInSurfaceRatios.mat'), 'infoPerSurfaceRatio', 'neighbours');

    end
end

