function unroll_OnlyApicalAndBasal(selpath)
%UNROLL_ONLYAPICALANDBASAL Summary of this function goes here
%   Detailed explanation goes here

    load(fullfile(selpath, '3d_layers_info.mat'), 'apicalLayer', 'basalLayer', 'glandOrientation', 'lumenImage', 'labelledImage');
    
    %% Creating image with a real size
    outputDirResults = strsplit(selpath, 'Results');
    zScaleFile = fullfile(outputDirResults{1}, 'Results', 'zScaleOfGland.mat');
    if exist(zScaleFile, 'file') > 0
        load(zScaleFile)
    else
        zScale = inputdlg('Insert z-scale of Gland');
        zScale = str2double(zScale{1});
        save(zScaleFile, 'zScale');
    end
    resizeImg = zScale;
    labelledImage_realSize = imresize3(labelledImage, resizeImg, 'nearest');
    insideGland = imresize3(double(labelledImage>0), resizeImg, 'nearest');
    insideGland = insideGland>0.75;
    labelledImage_realSize(insideGland == 0) = 0;
    
    save(fullfile(selpath, '3d_layers_info.mat'), 'labelledImage_realSize', '-append');
    
    
    %% -------------------------- APICAL -------------------------- %%
    disp('Apical');
    mkdir(fullfile(selpath, 'unrolledGlands', 'gland_SR_1'));
    [~, apicalAreaValidCells, apicalRotationsOriginal] = unrollTube(apicalLayer, fullfile(selpath, 'unrolledGlands', 'gland_SR_1'), fullfile(selpath, 'valid_cells.mat'), fullfile(selpath, '3d_layers_info.mat'));

    %% -------------------------- BASAL -------------------------- %%
    disp('Basal');
    mkdir(fullfile(selpath, 'unrolledGlands', 'gland_SR_basal'));
    unrollTube(basalLayer, fullfile(selpath, 'unrolledGlands', 'gland_SR_basal'), fullfile(selpath, 'valid_cells.mat'), fullfile(selpath, '3d_layers_info.mat'), apicalAreaValidCells, apicalRotationsOriginal);
end

