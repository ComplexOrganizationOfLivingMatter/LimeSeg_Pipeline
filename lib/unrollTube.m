function [samiraTable, areaOfValidCells, rotationsOriginal] = unrollTube(img3d_original,img3dComplete, outputDir, noValidCells, colours, rotationsOriginal)
%UNROLLTUBE Unroll Salivary Glands as tubes
%   The cylindrical layers (apical and basal) of the salivary gland were
%   mapped in the Cartesian plane for analysis using a cylindrical 
%   coordinates transformation.
    
    load(noValidCells);
    load(colours,'colours');
    colours = vertcat([1 1 1], colours);
       
    clearvars basalLayer apicalLayer labelledImage lumenImage
    
    closingPxAreas3D = 10;
    closingPxAreas2D = closingPxAreas3D;
    
    mkdir(outputDir)
    if exist(fullfile(outputDir, 'final3DImg.mat'), 'file')
        load(fullfile(outputDir, 'final3DImg.mat'));
    else
                
        %% Rotate the image to obtain all the glands with the same orientation
        [img3d_original] = rotateImg3(img3d_original,rotationsOriginal);
        [img3dComplete] = rotateImg3(img3dComplete, rotationsOriginal);
        
        [allX,allY,allZ]=ind2sub(size(img3dComplete),find(img3dComplete>0));
        img3d_originalCropped = img3d_original(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
        img3dComplete_Cropped = img3dComplete(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
        
        %% Check which side is the longest
        sizeImg3d = size(img3d_originalCropped);
        [~, indices] = sort(sizeImg3d);
        img3d_original = permute(img3d_originalCropped, indices);
        img3dComplete = permute(img3dComplete_Cropped, indices);
        tipValue = 21;
        img3d_original = addTipsImg3D(tipValue, img3d_original);
        img3dComplete = addTipsImg3D(tipValue, img3dComplete);
                
        clearvars img3d_originalCropped img3dComplete_Cropped
        
        img3d = double(img3d_original);
        img3dComplete = double(img3dComplete);
        validRegion = double(img3d_original>0);
        
        [validRegion] = imclose(validRegion>0, strel('sphere', closingPxAreas3D));
        
        % In e-cadhi salivary glands the lumen is too irregular, so we
        % can't use img3dComplete to empty spaces. We will use the image
        % itself.
        if contains(outputDir, 'E-cadh') && contains(outputDir, 'SR_1')
            img3d = fill0sWithCells(img3d .* double(validRegion), img3d, (img3dComplete>0 & validRegion)==0);
        else
            img3d = fill0sWithCells(img3d .* double(validRegion), img3dComplete, (img3dComplete>0 & validRegion)==0);
        end
        
        save(fullfile(outputDir, 'final3DImg.mat'), 'img3d', 'img3dComplete', '-v7.3');
    end

    %% Step 2: Get each 3D spherical line to a 2D line
    if exist(fullfile(outputDir, 'verticesInfo.mat'), 'file') == 0
        
        [cylindre2DImage, newVerticesNeighs2D, newVertices2D, centroids, ...
            validCellsFinal, borderCells, surfaceRatio, nameOfSimulation, ...
            areaOfValidCells, deployedImg, deployedImg3x, wholeImage, ...
            polygon_distribution, newNeighbours2D, newNeighbours2D_Checked] = mappCylindricalCoordinatesInto2D(img3d, img3dComplete, closingPxAreas2D, noValidCells, colours, outputDir);
        
    else
        load(fullfile(outputDir, 'allInfo.mat'));
        outputDirOri = outputDir;
        load(fullfile(outputDir, 'verticesInfo.mat'));
        outputDir = outputDirOri;
    end

    %% Export as samira table and other features
    if exist(fullfile(outputDir, 'samiraTable.mat'), 'file') == 0
        outputDirSplitted = strsplit(outputDir, 'Results');
        load(fullfile(outputDirSplitted{1}, 'Results', 'valid_cells.mat'))
        samiraTable = connectVerticesOf2D(deployedImg, newVerticesNeighs2D, newVertices2D, centroids, validCells, borderCells, surfaceRatio, outputDir, nameOfSimulation, deployedImg3x, img3d);
        save(fullfile(outputDir, 'samiraTable.mat'), 'samiraTable');
    else
        load(fullfile(outputDir, 'samiraTable.mat'));
    end
end

