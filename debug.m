files = dir('**/Salivary gland/**/Results/dividedGland/glandDividedInSurfaceRatios.mat');
for numFile = 1:length(files)
    load(fullfile(files(numFile).folder, 'glandDividedInSurfaceRatios.mat'))
    %figure; paint3D( ismember(imageOfSurfaceRatios{numPartition, 3}, validCells) .* imageOfSurfaceRatios{numPartition, 3}, [], colours);
    intermidScutoids = sum(infoPerSurfaceRatio{6, 4}.Scutoids) / size(infoPerSurfaceRatio{6, 4}, 1);
    finalScutoids = sum(infoPerSurfaceRatio{11, 4}.Scutoids) / size(infoPerSurfaceRatio{11, 4}, 1);
    logScutoids(numFile) = intermidScutoids/finalScutoids;
end