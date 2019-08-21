function [samiraTable, areaOfValidCells, rotationsOriginal] = unrollTube(img3d_original, outputDir, noValidCells, colours, apicalArea, apicalRotationsOriginal)
%UNROLLTUBE Unroll Salivary Glands as tubes
%   The cylindrical layers (apical and basal) of the salivary gland were
%   mapped in the Cartesian plane for analysis using a cylindrical 
%   coordinates transformation.
    
    load(noValidCells);
    load(colours);
    colours = vertcat([1 1 1], colours);
    
    %% Step 1: Creating image with its real size, in case it is necessary
    if exist('labelledImage_realSize', 'var')
        img3dComplete = labelledImage_realSize;
    else
        img3dComplete = labelledImage;
    end
    
    closingPxAreas3D = 10;
    closingPxAreas2D = closingPxAreas3D;
    
    mkdir(outputDir)
    if exist(fullfile(outputDir, 'final3DImg.mat'), 'file')
        load(fullfile(outputDir, 'final3DImg.mat'));
        if exist('apicalRotationsOriginal', 'var') ~= 0
            rotationsOriginal = apicalRotationsOriginal;
        else
            load(fullfile(outputDir, 'apicalRotationsOriginal.mat'));
        end
    else
        outputDirResults = strsplit(outputDir, 'Results');
        zScaleFile = fullfile(outputDirResults{1}, 'Results', 'zScaleOfGland.mat');
        if exist(zScaleFile, 'file') > 0
            load(zScaleFile)
        else
            zScale = inputdlg('Insert z-scale of Gland');
            zScale = str2double(zScale{1});
            save(zScaleFile, 'zScale');
        end
        
        if exist('labelledImage_realSize', 'var')
            resizeImg = 1;
        else
            resizeImg = 1/zScale;
        end
        
        %% Rotate the image to obtain all the glands with the same orientation
        if exist('apicalRotationsOriginal', 'var') == 0
            [img3d_original, rotationsOriginal] = rotateImg3(img3d_original);
            [img3dComplete] = rotateImg3(img3dComplete, rotationsOriginal);
            save(fullfile(outputDir, 'apicalRotationsOriginal.mat'), 'rotationsOriginal');
        else
            [img3d_original] = rotateImg3(img3d_original, apicalRotationsOriginal);
            [img3dComplete] = rotateImg3(img3dComplete, apicalRotationsOriginal);
            rotationsOriginal = apicalRotationsOriginal;
        end
        
        [allX,allY,allZ]=ind2sub(size(img3dComplete),find(img3dComplete>0));
        img3d_originalCropped = img3d_original(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
        img3dComplete_Cropped = img3dComplete(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
        
        % Check which side is the longest
        %img3d_closed = fill0sWithCells(img3d_original, img3dComplete, (img3dComplete & imclose(img3d_original>0, strel('sphere', closingPxAreas))) == 0);
        %neighbours = getNeighboursFromFourProjectedPlanesFrom3Dgland(img3d_closed, colours);
        %neighbours = checkPairPointCloudDistanceCurateNeighbours(img3d_closed, neighbours);

        sizeImg3d = size(img3d_originalCropped);
        [~, indices] = sort(sizeImg3d);
        img3d_original = permute(img3d_originalCropped, indices);
        img3dComplete = permute(img3dComplete_Cropped, indices);
        tipValue = 21;
        img3d_original = addTipsImg3D(tipValue, img3d_original);
        img3dComplete = addTipsImg3D(tipValue, img3dComplete);
        %[verticesInfo] = getVertices3D(img3d_original, neighbours);
        %vertices3D = vertcat(verticesInfo.verticesPerCell{:});

%         [colours] = exportAsImageSequence_WithVertices(img3d, outputDir, colours, -1, vertices3D);
%         figure; paint3D(img3d_original, [], colours);
%         for numVertex = 1:size(vertices3D, 1)
%             hold on; plot3(vertices3D(numVertex, 1), vertices3D(numVertex, 2), vertices3D(numVertex, 3), 'rx')
%         end

        %vertices3D_Neighbours = verticesInfo.verticesConnectCells;
        %vertices3D_Neighbours(cellfun(@isempty, verticesInfo.verticesPerCell), :) = [];
        %cellNumNeighbours = cellfun(@length, neighbours);
        
        img3d = double(img3d_original);
        img3dComplete = double(img3dComplete);
        imgSize = round(size(img3d_original)/resizeImg);
        img3d = double(imresize3(img3d_original, imgSize, 'nearest'));
        
        img3dComplete = double(imresize3(img3dComplete, imgSize, 'nearest'));
        
        validRegion = double(imresize3(img3d_original, imgSize)>0);
        
        [validRegion] = imclose(validRegion>0, strel('sphere', closingPxAreas3D));
                    
        img3d = fill0sWithCells(img3d .* double(validRegion), img3dComplete, (img3dComplete>0 & validRegion)==0);
        %vertices3D = round(vertices3D / resizeImg);
        mkdir(outputDir);
        %save(fullfile(outputDir, 'final3DImg.mat'), 'img3d', 'img3dComplete', 'vertices3D_Neighbours', 'vertices3D', 'cellNumNeighbours', 'neighbours', '-v7.3');
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

