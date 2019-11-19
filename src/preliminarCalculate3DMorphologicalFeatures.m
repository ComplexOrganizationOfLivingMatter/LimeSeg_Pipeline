function [] = preliminarCalculate3DMorphologicalFeatures()
%CALCULATE3DMORPHOLOGICALFEATURES Summary of this function goes here
addpath(genpath('lib'))
addpath(genpath('src'))

clear all
%   Detailed explanation goes here
folderName = uigetdir('data');
if contains(folderName, '\Salivary gland\')
    files = dir(fullfile(folderName, '**/Results/3d_layers_info.mat'));
else
    files = dir(fullfile('**/Salivary gland/', folderName, '/**/Results/3d_layers_info.mat'));
end
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'incomplete') == 0, {files.folder});
files = files(nonDiscardedFiles);

totalMeanFeatures = [];
totalStdFeatures = [];
allGlands = [];
allLumens = [];
allFilesName = [];

for numFiles=1:length(files)
%     if exist(fullfile(files(numFiles).folder, 'preliminarMorphological3dFeatures.mat'), 'file') == 0
        files(numFiles).folder
        load(fullfile(files(numFiles).folder, '3d_layers_info.mat'), 'lumenImage', 'labelledImage');
            %% Creating image with a real size
            outputDirResults = strsplit(files(numFiles).folder, 'Results');
        zScaleFile = fullfile(files(numFiles).folder, 'zScaleOfGland.mat');
            if exist(zScaleFile, 'file') > 0
                load(zScaleFile)
            else
                zScale = inputdlg('Insert z-scale of Gland');
                zScale = str2double(zScale{1});
                save(zScaleFile, 'zScale');
            end
            resizeImg = zScale;

            labelledImage_realSize = imresize3(labelledImage, resizeImg, 'nearest');
            lumenImage_realSize = imresize3(double(lumenImage), resizeImg, 'nearest');        
            load(fullfile(files(numFiles).folder, 'valid_cells.mat'), 'validCells');
            
            %% Obtain layers on its real 3D size
        basalLayer = getBasalFrom3DImage(labelledImage_realSize, lumenImage_realSize, 0, labelledImage_realSize == 0 & lumenImage_realSize == 0);
        [apicalLayer] = getApicalFrom3DImage(lumenImage_realSize, labelledImage_realSize);
        %% Basal features
        basal_area_3dcells=cell2mat(struct2cell(regionprops(basalLayer,'Area'))).';
        basal_area_3dcells = basal_area_3dcells(validCells);
        
        %% Apical features
        apical_area_3dcells=cell2mat(struct2cell(regionprops(apicalLayer,'Area'))).';
        apical_area_3dcells = apical_area_3dcells(validCells);
        
        %% Extract each cell and calculate 3D features
        [cells3dFeatures] = extract3dDescriptors(labelledImage_realSize, validCells');
        cells3dFeatures = horzcat(cells3dFeatures, table(apical_area_3dcells, basal_area_3dcells));

        %% Lumen features
        [lumen3dFeatures] = extract3dDescriptors(lumenImage_realSize, 1);
        lumen3dFeatures.ID_Cell = 'Lumen';
        lumen3dFeatures.basal_area_cells = -1;
        lumen3dFeatures.apical_area_cells = -1;
        gland3dFeatures.percentageScutoids = -1;
        gland3dFeatures.totalNeighs = -1;

        %% Global Gland
        % We need calculate thickness of the glands or number of cell in
        % transversal axis
        [gland3dFeatures] = extract3dDescriptors(labelledImage_realSize>0, 1);

        gland3dFeatures.ID_Cell = 'Gland';
        gland3dFeatures.basal_area_cells = -1;
        gland3dFeatures.apical_area_cells = -1;
        gland3dFeatures.percentageScutoids = -1;
        gland3dFeatures.totalNeighs = -1;
        
        %% Save variables and export to excel
        writetable(cells3dFeatures,fullfile(files(numFiles).folder,'3dFeatures_LimeSeg3DSegmentation.xls'), 'Range','B2');
        save(fullfile(files(numFiles).folder, 'preliminarMorphological3dFeatures.mat'), 'cells3dFeatures', 'gland3dFeatures', 'lumen3dFeatures');
%     else
%         load(fullfile(files(numFiles).folder, 'preliminarMorphological3dFeatures.mat'));
%     end
    
    %% Calculate mean and std of 3D features
    meanFeatures = varfun(@(x) mean(x),cells3dFeatures(:, 2:end));
    Fun_numberCells =  sum(height(cells3dFeatures(:, 1)));
    stdFeatures = varfun(@(x) std(x),cells3dFeatures(:, 2:end));
    Fun_SurfaceRatio = meanFeatures.Fun_basal_area_3dcells/meanFeatures.Fun_apical_area_3dcells;
    meanFeatures = horzcat(meanFeatures, table(Fun_SurfaceRatio, Fun_numberCells));
    
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

    save(fullfile(folderName, 'global_3dFeatures.mat'), 'totalMeanFeatures','totalStdFeatures', 'allLumens', 'allGlands')
    writetable([allFilesName,totalMeanFeatures,totalStdFeatures, allGlands, allLumens], fullfile(folderName,'global_3dFeatures.xls'),'Range','B2');
end

