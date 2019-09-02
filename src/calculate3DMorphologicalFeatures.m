function [] = calculate3DMorphologicalFeatures(folderName)
%CALCULATE3DMORPHOLOGICALFEATURES Summary of this function goes here
%   Detailed explanation goes here

files = dir(fullfile('**/data/Salivary gland/', folderName, '/**/Results/3d_layers_info.mat'));
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
files = files(nonDiscardedFiles);

selpath = dir(fullfile('**/data/Salivary gland/', folderName));

totalMeanFeatures = [];
totalStdFeatures = [];
allGlands = [];
allLumens = [];
allFilesName = [];

for numFiles=1:length(files)
    files(numFiles).folder
    if exist(fullfile(files(numFiles).folder, 'morphological3dFeatures.mat'), 'file') == 0
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
        save(fullfile(files(numFiles).folder, 'morphological3dFeatures.mat'), 'cells3dFeatures', 'gland3dFeatures', 'lumen3dFeatures');
    else
        load(fullfile(files(numFiles).folder, 'morphological3dFeatures.mat'));
    end
    
    %% Calculate mean and std of 3D features
    meanFeatures = varfun(@(x) mean(x),cells3dFeatures(:, 2:end));
    stdFeatures = varfun(@(x) std(x),cells3dFeatures(:, 2:end));
    
    totalMeanFeatures = vertcat(totalMeanFeatures, meanFeatures);
    totalStdFeatures = vertcat(totalStdFeatures, stdFeatures);
    allGlands = vertcat(allGlands, gland3dFeatures);
    allLumens = vertcat(allLumens, lumen3dFeatures);
    
    fileName = strsplit(files(numFiles).folder, {'/','\'});
    fileName = convertCharsToStrings(strjoin({fileName{1,end-2},fileName{1,end-1}}, ' '));
    allFilesName = [allFilesName ; fileName];
end


allFilesName = table(allFilesName, 'VariableNames', {'ID_Glands'});

allGlands.Properties.VariableNames = cellfun(@(x) strcat('Gland_', x), allGlands.Properties.VariableNames, 'UniformOutput', false);
allLumens.Properties.VariableNames = cellfun(@(x) strcat('Lumen_', x), allLumens.Properties.VariableNames, 'UniformOutput', false);
totalMeanFeatures.Properties.VariableNames = cellfun(@(x) strcat('AverageCell_', x(5:end)), totalMeanFeatures.Properties.VariableNames, 'UniformOutput', false);
totalStdFeatures.Properties.VariableNames = cellfun(@(x) strcat('STDCell_', x(5:end)), totalStdFeatures.Properties.VariableNames, 'UniformOutput', false);

save(fullfile(selpath(1).folder, 'global_3dFeatures.mat'), 'totalMeanFeatures','totalStdFeatures', 'allLumens', 'allGlands')
writetable([allFilesName,totalMeanFeatures,totalStdFeatures, allGlands, allLumens], fullfile(selpath(1).folder,'global_3dFeatures.xls'),'Range','B2');
end
