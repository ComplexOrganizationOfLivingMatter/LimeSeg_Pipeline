function unroll_OnlyApicalAndBasal(selpath, testing)
%UNROLL_ONLYAPICALANDBASAL Perform unrolling only to the apical and basal
%surface of a S. Gland
%   
    if exist(fullfile(selpath, '/unrolledGlands/gland_SR_1/verticesInfo.mat'), 'file')==0 || exist(fullfile(selpath, '/unrolledGlands/gland_SR_basal/verticesInfo.mat'), 'file')==0
        variablesOfFile = who('-file', fullfile(selpath, '3d_layers_info.mat'));

        if exist(fullfile(selpath, 'realSize3dLayers.mat'),'file')
            load(fullfile(selpath, 'realSize3dLayers.mat'), 'lumenImage_realSize', 'labelledImage_realSize');
        else
            if any(cellfun(@(x) isequal('labelledImage_realSize', x), variablesOfFile)) == 0 || ~contains(lower(selpath), 'flatten')
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

                if contains(lower(selpath), 'flatten')
                    [labelledImage] = flattenMutantGland(apicalLayer, basalLayer, labelledImage, lumenImage);
                end

                labelledImage_realSize = imresize3(labelledImage, resizeImg, 'nearest');
                lumenImage_realSize = imresize3(double(lumenImage), resizeImg, 'nearest');
            %     insideGland = imresize3(double(labelledImage>0), resizeImg, 'nearest');
            %     insideGland = insideGland>0.75;
            %     labelledImage_realSize(insideGland == 0) = 0;

                save(fullfile(selpath, '3d_layers_info.mat'), 'labelledImage_realSize', 'lumenImage_realSize', '-append');
                [msgstr] = lastwarn;
                if contains(msgstr, 'was not saved')
                    ME = MException('MyComponent:noVariableSaved', ...
                    'A variable was not saved. Please use -v7.3');
                    throw(ME)
                end
            else 
                load(fullfile(selpath, '3d_layers_info.mat'), 'lumenImage_realSize', 'labelledImage_realSize');
            end
        end

        labelledImage_realSize=uint16(labelledImage_realSize);
        lumenImage_realSize=uint16(lumenImage_realSize);
        %% Obtain layers on its real 3D size
        basalLayer = getBasalFrom3DImage(labelledImage_realSize, lumenImage_realSize, 0, labelledImage_realSize == 0 & lumenImage_realSize == 0);
        [apicalLayer] = getApicalFrom3DImage(lumenImage_realSize, labelledImage_realSize);

        if contains(lower(selpath), 'e-cadhi')
           variablesOfFile = who('-file', fullfile(selpath, 'realSize3dLayers.mat'));
           if any(cellfun(@(x) isequal('labelledImage_realSizeFlatten', x), variablesOfFile))
               load(fullfile(selpath, 'realSize3dLayers.mat'), 'labelledImage_realSizeFlatten');
           else
               [labelledImage_realSizeFlatten] = flattenMutantGlandUsingLateralCellWalls(labelledImage_realSize, lumenImage_realSize);
               emptyVoxels2Fill = labelledImage_realSizeFlatten(apicalLayer>0)==0;
               labelledImage_realSizeFlatten(emptyVoxels2Fill)=apicalLayer(emptyVoxels2Fill);
               save(fullfile(selpath, 'realSize3dLayers.mat'), 'labelledImage_realSizeFlatten','-append');
           end
           basalLayer = getBasalFrom3DImage(labelledImage_realSizeFlatten, lumenImage_realSize, 0, labelledImage_realSizeFlatten == 0 & lumenImage_realSize == 0);
        end

    % 
    %     figure; paint3D(apicalLayer)
    %     figure; paint3D(basalLayer)
   
        %% -------------------------- APICAL UNROLL-------------------------- %%
        disp('Apical');
        mkdir(fullfile(selpath, 'unrolledGlands', 'gland_SR_1'));
        %% Get gland rotation using the gland lumen
        if exist(fullfile(selpath, 'unrolledGlands', 'gland_SR_1','apicalRotationsOriginal.mat'), 'file') == 0
            [rotatedLumen, rotationsOriginal]= rotateImg3(lumenImage_realSize);
            save(fullfile(selpath, 'unrolledGlands', 'gland_SR_1','apicalRotationsOriginal.mat'),'rotationsOriginal');
        else
            load(fullfile(selpath, 'unrolledGlands', 'gland_SR_1','apicalRotationsOriginal.mat'),'rotationsOriginal');
        end
        unrollTube(apicalLayer,labelledImage_realSizeFlatten, fullfile(selpath, 'unrolledGlands', 'gland_SR_1'), fullfile(selpath, 'valid_cells.mat'), fullfile(selpath, '3d_layers_info.mat'),rotationsOriginal);

        %% -------------------------- BASAL UNROLL-------------------------- %%
        disp('Basal');
        mkdir(fullfile(selpath, 'unrolledGlands', 'gland_SR_basal'));
        unrollTube(basalLayer,labelledImage_realSizeFlatten, fullfile(selpath, 'unrolledGlands', 'gland_SR_basal'), fullfile(selpath, 'valid_cells.mat'), fullfile(selpath, '3d_layers_info.mat'), rotationsOriginal);
    end
end