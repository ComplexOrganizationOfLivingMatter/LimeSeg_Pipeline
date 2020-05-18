function [basalLayer,apicalLayer,labelledImage_realSize,lumenImage_realSize]=resizeTissue(numFile,files)

if exist(fullfile(files(numFile).folder, 'realSize3dLayers.mat'), 'file') == 0
    
    %% Step 1: Creating image with its real size
    load(fullfile(files(numFile).folder, 'zScaleOfGland'))
    load(fullfile(files(numFile).folder, '3d_layers_info.mat'))%, 'labelledImage_realSize', 'lumenImage_realSize');
    
    if contains(lower(files(numFile).folder), 'flatten')
        [labelledImage] = flattenMutantGland(apicalLayer, basalLayer, labelledImage, lumenImage);
    end
    
    if contains(lower(files(numFile).folder), 'oldmethod')
        labelledImage_realSize = imresize3(labelledImage, zScale, 'nearest');
        lumenImage_realSize = imresize3(double(lumenImage), zScale, 'nearest');
    else
        labelledImage_realSize  = imresize3(labelledImage, [1024 1024 zScale*92], 'nearest');
        lumenImage_realSize  = imresize3(double(lumenImage), [1024 1024 zScale*92], 'nearest');
        lumenImage_realSize  = logical(lumenImage_realSize );    
    end
    
    basalLayer = getBasalFrom3DImage(labelledImage_realSize, lumenImage_realSize, 0, labelledImage_realSize == 0 & lumenImage_realSize == 0);
    [apicalLayer] = getApicalFrom3DImage(lumenImage_realSize, labelledImage_realSize);
    
    layers3d=[{apicalLayer},{basalLayer}];
    
    %% Step 2: Getting the perimeter of the basal and apical layers of the image with its real size
    for iteration=1:2
        
        img3dComplete = labelledImage_realSize;
        img3d_original = cell2mat(layers3d(iteration));
        
        [allX,allY,allZ]=ind2sub(size(img3dComplete),find(img3dComplete>0));
        img3d_originalCropped = img3d_original(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
        img3dComplete_Cropped = img3dComplete(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
        
        sizeImg3d = size(img3d_originalCropped);
        [~, indices] = sort(sizeImg3d);
        img3d_original = permute(img3d_originalCropped, indices);
        img3dComplete = permute(img3dComplete_Cropped, indices);
        tipValue = 20;
        img3d_original = addTipsImg3D(tipValue, img3d_original);
        img3dComplete = addTipsImg3D(tipValue, img3dComplete);
        
        img3dComplete = double(img3dComplete);
        img3d = double(img3d_original);
        
        if iteration == 1
            [apicalLayer] = calculatePerimOf3DImage(img3d, img3dComplete);
        else
            [basalLayer] = calculatePerimOf3DImage(img3d, img3dComplete);
        end
        
    end
    save(fullfile(files(numFile).folder, 'realSize3dLayers.mat'), 'labelledImage_realSize','lumenImage_realSize','apicalLayer','basalLayer', '-v7.3');
    
else
    load(fullfile(files(numFile).folder, 'realSize3dLayers.mat'))
end

end