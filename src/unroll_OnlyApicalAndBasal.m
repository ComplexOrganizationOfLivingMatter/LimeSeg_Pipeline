function unroll_OnlyApicalAndBasal(selpath, testing)
%UNROLL_ONLYAPICALANDBASAL Perform unrolling only to the apical and basal
%surface of a S. Gland
%   

    variablesOfFile = who('-file', fullfile(selpath, '3d_layers_info.mat'));
    
    if any(cellfun(@(x) isequal('labelledImage_realSize', x), variablesOfFile))
        load(fullfile(selpath, '3d_layers_info.mat'), 'lumenImage_realSize', 'labelledImage_realSize');
    else
        load(fullfile(selpath, '3d_layers_info.mat'), 'lumenImage', 'labelledImage', 'apicalLayer', 'basalLayer');
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
        
            
        if contains(lower(selpath), 'e-cadh')
            [labelledImage] = flattenMutantGland(apicalLayer, basalLayer, labelledImage);
        end
        
        labelledImage_realSize = imresize3(labelledImage, resizeImg, 'nearest');
        lumenImage_realSize = imresize3(double(lumenImage), resizeImg, 'nearest');
    %     insideGland = imresize3(double(labelledImage>0), resizeImg, 'nearest');
    %     insideGland = insideGland>0.75;
    %     labelledImage_realSize(insideGland == 0) = 0;
        
        save(fullfile(selpath, '3d_layers_info.mat'), 'labelledImage_realSize', 'lumenImage_realSize', '-append');
    end
    
    %% Obtain layers on its real 3D size
    basalLayer = getBasalFrom3DImage(labelledImage_realSize, lumenImage_realSize, 0, labelledImage_realSize == 0 & lumenImage_realSize == 0);
    [apicalLayer] = getApicalFrom3DImage(lumenImage_realSize, labelledImage_realSize);

% 
%     figure; paint3D(apicalLayer)
%     figure; paint3D(basalLayer)
    
    %% -------------------------- APICAL -------------------------- %%
    disp('Apical');
    mkdir(fullfile(selpath, 'unrolledGlands', 'gland_SR_1'));
    [~, apicalAreaValidCells, apicalRotationsOriginal] = unrollTube(apicalLayer, fullfile(selpath, 'unrolledGlands', 'gland_SR_1'), fullfile(selpath, 'valid_cells.mat'), fullfile(selpath, '3d_layers_info.mat'));

    %% -------------------------- BASAL -------------------------- %%
    disp('Basal');
    mkdir(fullfile(selpath, 'unrolledGlands', 'gland_SR_basal'));
    unrollTube(basalLayer, fullfile(selpath, 'unrolledGlands', 'gland_SR_basal'), fullfile(selpath, 'valid_cells.mat'), fullfile(selpath, '3d_layers_info.mat'), apicalAreaValidCells, apicalRotationsOriginal);
end

