function [] = calculate3DMorphologicalFeatures(folderName)
%CALCULATE3DMORPHOLOGICALFEATURES Summary of this function goes here
%   Detailed explanation goes here

files = dir(fullfile('**/data/Salivary gland/', folderName, '/**/Results/3d_layers_info.mat'));
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
files = files(nonDiscardedFiles);

totalMeanFeatures = [];
totalStdFeatures = [];
allGlands = [];
allLumens = [];
allFilesName = [];

for numFiles=1:length(files)
    load(fullfile(files(numFiles).folder, '3d_layers_info.mat'), 'labelledImage_realSize', 'lumenImage_realSize');
    load(fullfile(files(numFiles).folder, 'valid_cells.mat'), 'validCells');
    
    %% Extract each cell and calculate 3D features
    [cells3dFeatures] = extract3dDescriptors(labelledImage_realSize, validCells');
    
    %% Lumen features
    [lumen3dFeatures] = extract3dDescriptors(lumenImage_realSize, 1);
    lumen3dFeatures.ID_Cell = 'Lumen';
    
    %% Global Gland
    % We need calculate thickness of the glands or number of cell in
    % transversal axis
    [gland3dFeatures] = extract3dDescriptors(labelledImage_realSize>0, 1);
    gland3dFeatures.ID_Cell = 'Gland';
    
    allFeatures = vertcat(cells3dFeatures, gland3dFeatures, lumen3dFeatures);
    %% Save variables and export to excel
    writetable(allFeatures,fullfile(files(numFiles).folder,'3dFeatures_LimeSeg3DSegmentation.xls'), 'Range','B2');
    save(fullfile(files(numFiles).folder, 'morphological3dFeatures.mat'), 'allFeatures')
    
    %% Calculate mean and std of 3D features
    meanFeatures = varfun(@(x) mean(x),cells3dFeatures(2:end, :));
    stdFeatures = varfun(@(x) std(x),cells3dFeatures(2:end, :));
    
    totalMeanFeatures(end+1) = meanFeatures;
    totalStdFeatures(end+1) = stdFeatures;
    allGlands(end+1) = gland3dFeatures;
    allLumens(end+1) = lumen3dFeatures;
    
    fileName = strsplit(files(numFiles).folder, {'/','\'});
    fileName = convertCharsToStrings(strjoin({fileName{1,end-2},fileName{1,end-1}}, ' '));
    allFilesName = [allFilesName ; fileName];
end

selpath = dir(fullfile('**/data/Salivary gland/', folderName));

allFilesName = table(allFilesName, 'VariableNames', {'ID_Glands'});
totalStdFeatures.Properties.VariableNames = {'stdVolume','stdEquivDiameter','stdPrincipalAxisLength','stdConvexVolume', 'stdSolidity', 'stdSurfaceArea', 'stdAspectRatio', 'stdSphericity', 'stdNormalizedVolume'};

save(fullfile(selpath(1).folder, 'global_3dFeatures.mat'), 'totalMeanFeatures','totalStdFeatures')
writetable([allFilesName,totalMeanFeatures,totalStdFeatures, allGlands, allLumens], fullfile(selpath(1).folder,'global_3dFeatures.xls'),'Range','B2');
end

