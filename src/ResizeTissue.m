function [basalLayer,apicalLayer,labelledImage_realSize,lumenImage_realSize]=ResizeTissue(numFile,files,labelledImage,lumenImage,zScaleGland)

%% Step 1: Creating image with its real size, in case it is necessary
load(fullfile(files(numFile).folder, 'zScaleOfGland'))
if exist(fullfile(files(numFile).folder, 'realSize3dLayers.mat'), 'file') == 0
    labelledImage_realSize = imresize3(labelledImage, zScale, 'nearest');
    lumenImage_realSize = imresize3(double(lumenImage), zScale, 'nearest');
    
    basalLayer = getBasalFrom3DImage(labelledImage_realSize, lumenImage_realSize, 0, labelledImage_realSize == 0 & lumenImage_realSize == 0);
    [apicalLayer] = getApicalFrom3DImage(lumenImage_realSize, labelledImage_realSize);
    save(fullfile(files(numFile).folder, 'realSize3dLayers.mat'), 'labelledImage_realSize','lumenImage_realSize','apicalLayer','basalLayer', '-v7.3');
else
    load(fullfile(files(numFile).folder, 'realSize3dLayers.mat'))
end


layers=[{apicalLayer},{basalLayer}];

for iteration=1:2
    if exist(fullfile(files(numFile).folder, 'final3DImg.mat'), 'file')
        load(fullfile(files(numFile).folder, 'final3DImg.mat'));
    else
        
        if exist('labelledImage_realSize', 'var')
            resizeImg = 1;
        else
            resizeImg = 1/zScale;
        end
        
        img3dComplete = labelledImage_realSize;
        img3d_original = cell2mat(layers(iteration));
        
                %% Rotate the image to obtain all the glands with the same orientation
        if iteration == 1
            [img3d_original, rotationsOriginal] = rotateImg3(img3d_original);
            [img3dComplete] = rotateImg3(img3dComplete, rotationsOriginal);
            save(fullfile(files(numFile).folder, 'apicalRotationsOriginal.mat'), 'rotationsOriginal');
        else
            load(fullfile(files(numFile).folder, 'apicalRotationsOriginal.mat'), 'rotationsOriginal');
            [img3d_original] = rotateImg3(img3d_original, rotationsOriginal);
            [img3dComplete] = rotateImg3(img3dComplete, rotationsOriginal);

        end
        
        
        [allX,allY,allZ]=ind2sub(size(img3dComplete),find(img3dComplete>0));
        img3d_originalCropped = img3d_original(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
        img3dComplete_Cropped = img3dComplete(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
        
        sizeImg3d = size(img3d_originalCropped);
        [~, indices] = sort(sizeImg3d);
        img3d_original = permute(img3d_originalCropped, indices);
        img3dComplete = permute(img3dComplete_Cropped, indices);
        tipValue = 21;
        img3d_original = addTipsImg3D(tipValue, img3d_original);
        img3dComplete = addTipsImg3D(tipValue, img3dComplete);
        
        img3dComplete = double(img3dComplete);
        imgSize = round(size(img3d_original)/resizeImg);
        img3d = double(imresize3(img3d_original, imgSize, 'nearest'));
        
        img3dComplete = double(imresize3(img3dComplete, imgSize, 'nearest'));
        
        
        closingPxAreas3D = 10;
        
        validRegion = double(imresize3(img3d_original, imgSize)>0);
        [validRegion] = imclose(validRegion>0, strel('sphere', closingPxAreas3D));
        img3d = fill0sWithCells(img3d .* double(validRegion), img3dComplete, (img3dComplete & validRegion)==0);
        
        if iteration == 1
            [apicalLayer] = calculatePerimOf3DImage(img3d, img3dComplete);
        else
            [basalLayer] = calculatePerimOf3DImage(img3d, img3dComplete);
        end
    end
end
    
    save(fullfile(files(numFile).folder, 'final3DImg.mat'), 'apicalLayer','basalLayer');

end