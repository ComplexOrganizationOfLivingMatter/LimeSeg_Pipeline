function [cylindre2DImage, newVerticesNeighs2D, newVertices2D, centroids, validCellsFinal, borderCells, surfaceRatio, nameOfSimulation, areaOfValidCells, deployedImg, deployedImg3x, wholeImage, polygon_distribution, newNeighbours2D, newNeighbours2D_Checked] = mappCylindricalCoordinatesInto2D(img3d, img3dComplete, closingPxAreas2D, noValidCells, colours, outputDir)
%MAPPCYLINDRICALCOORDINATESINTO2D Summary of this function goes here
%   Detailed explanation goes here

    pixelSizeThreshold = 10;

    imgFinalCoordinates=cell(size(img3d,3),1);
    imgFinalCoordinates3x=cell(size(img3d,3),1);
    for coordZ = 1 : size(img3d,3)
        if sum(sum(img3d(:, :, coordZ) > 0)) < pixelSizeThreshold || sum(sum(img3d(:, :, coordZ))) < pixelSizeThreshold
            continue
        end
        
        closedZFrame = imclose(img3d(:, :, coordZ)>0, strel('disk', round(closingPxAreas2D)));
        img3d(:, :, coordZ) = fill0sWithCells(img3d(:, :, coordZ), img3dComplete(:, :, coordZ), closedZFrame==0);

        %% Remove pixels surrounding the boundary
        [filledImage] = createCompleteSection(img3d, coordZ, img3dComplete);

        %% Create perim
        [orderedLabels] = perim2line(filledImage, img3d, img3dComplete, coordZ);

        imgFinalCoordinates3x{coordZ} = repmat(orderedLabels, 1, 3);

        %If we find that a cell does not continue in the other side and
        %finish on one side, we put the same pixel on the other side.
        if orderedLabels(1) ~= orderedLabels(end)
            orderedLabels(end+1) = orderedLabels(1);
        end
        %length(orderedLabels)
        imgFinalCoordinates{coordZ} = orderedLabels;
    end

    %% Reconstruct deployed img
    ySize=max(cellfun(@length, imgFinalCoordinates3x));
    deployedImg3x = zeros(size(img3d,3),ySize);
    deployedImg = zeros(size(img3d,3),ySize);

    for coordZ = 1 : size(img3d,3)
        rowOfCoord3x = imgFinalCoordinates3x{coordZ};
        rowOfCoord = imgFinalCoordinates{coordZ};

        nEmptyPixels3x = 0;
        if length(rowOfCoord3x) < ySize
            nEmptyPixels3x = floor((ySize - length(rowOfCoord3x)) / 2);
            nEmptyPixels = floor((ySize - length(rowOfCoord)) / 2);
        end
        deployedImg3x(coordZ, 1 + nEmptyPixels3x : length(rowOfCoord3x) + nEmptyPixels3x) = rowOfCoord3x;
        deployedImg(coordZ, 1 + nEmptyPixels : length(rowOfCoord) + nEmptyPixels) = rowOfCoord;
    end

    %% Getting correct border cells, valid cells and no valid cells
    %cylindre2DImage = fillEmptySpacesByWatershed2D(deployedImg, imclose(deployedImg>0, strel('disk', 20)) == 0 , colours);
    validCellsFinal  = setdiff(1:max(deployedImg(:)), noValidCells);
    deployedImg = fill0sWithCells(deployedImg, deployedImg, imfill(ismember(deployedImg, validCellsFinal)>0, 'holes')==0);
    cylindre2DImage = deployedImg;
    %          figure;imshow(deployedImg,colours)
    deployedImg3x = fill0sWithCells(deployedImg3x, deployedImg3x, imfill(ismember(deployedImg3x, validCellsFinal)>0, 'holes')==0);
    [wholeImage] = fillEmptySpacesByWatershed2D(deployedImg3x, imclose(deployedImg3x>0, strel('disk', 3)) == 0 , colours);

    %     figure;imshow(finalImage,colours)

    %% We only keep the cells in the middle
    relabelFinalImage = bwlabel(wholeImage,4);
    labelsFinal = unique(relabelFinalImage(deployedImg>0));
    midSectionImage = wholeImage;
    midSectionImage(~ismember(relabelFinalImage,labelsFinal))=0;

    %% We keep the valid cells from that middle image
    finalImageWithValidCells = ismember(cylindre2DImage, validCellsFinal).*cylindre2DImage;
    %     figure;imshow(finalImageWithValidCells,colours)

    if exist('outputDir', 'var')
        h = figure ('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
        imshow(midSectionImage+1, colours);
        hold on;
        midSectionNewLabels = bwlabel(midSectionImage, 4);
        centroids = regionprops(midSectionNewLabels, 'Centroid');
        centroids = round(vertcat(centroids.Centroid));

        ax = get(h, 'Children');
        set(ax,'Units','normalized')
        set(ax,'Position',[0 0 1 1])
        for numCentroid = 1:size(centroids, 1)
            labelSeed = midSectionImage(midSectionNewLabels == numCentroid);
            labelSeed = labelSeed(1);
            if mean(colours(labelSeed, :)) < 0.4
                text(ax, centroids(numCentroid, 1), centroids(numCentroid, 2), num2str(labelSeed), 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 6);
            else
                text(ax, centroids(numCentroid, 1), centroids(numCentroid, 2), num2str(labelSeed), 'HorizontalAlignment', 'center', 'FontSize', 6);
            end
        end
        h.InvertHardcopy = 'off';
        saveas(h, fullfile(outputDir, 'img_MidSection.tif'));
        imwrite(finalImageWithValidCells+1, colours, fullfile(outputDir, 'img_MidSection_ValidCells.tif'));
        imwrite(wholeImage+1, colours, fullfile(outputDir, 'img_WholeImage.tif'));
        imwrite(deployedImg+1, colours, fullfile(outputDir, 'img_original.tif'));
    else
        centroids = [];
    end

    %% Calculating surface ratio
    areaValidCellsImg = deployedImg .* ismember(deployedImg, validCellsFinal);
    areaOfValidCells = sum(areaValidCellsImg(:) > 0);

    if exist('apicalArea', 'var') == 0
        surfaceRatio = 1;
    else
        surfaceRatio = areaOfValidCells / apicalArea;
    end

    %% Connect vertices to obtain an image from the vertices
    newNeighbours2D = calculateNeighbours(deployedImg);
    newNeighbours2D_Checked = checkPairPointCloudDistanceCurateNeighbours(img3d, newNeighbours2D);

    newVertices2D = getVertices(deployedImg, newNeighbours2D_Checked);
    newVerticesNeighs2D = vertcat(newVertices2D.verticesConnectCells);
    newVerticesNeighs2D_empty = cellfun(@isempty, newVertices2D.verticesPerCell);
    newVerticesNeighs2D(newVerticesNeighs2D_empty, :) = [];
    newVertices2D = vertcat(newVertices2D.verticesPerCell{:});

    toGetBorderCells = cylindre2DImage+1;
    background = imopen(cylindre2DImage==0, strel('sphere', 1));
    toGetBorderCells(~background & cylindre2DImage == 0) = 0;
    backgroundNeighs = calculateNeighbours(toGetBorderCells);
    borderCells = backgroundNeighs{1} - 1;
    borderCells = intersect(validCellsFinal, borderCells);
    if exist('apicalArea', 'var') == 0
        nameOfSimulation = 'Apical';
    else
        nameOfSimulation = 'Basal';
    end

    occurrencesOfCells = tabulate(newVerticesNeighs2D(:));
    side_cells = occurrencesOfCells(:, 2);
    [polygon_distribution] = calculate_polygon_distribution(side_cells, validCellsFinal);
    polygon_distribution_T = cell2table(polygon_distribution(2:end, :));
    polygon_distribution_T.Properties.VariableNames = polygon_distribution(1, :);
    
    if exist('outputDir', 'var')
        writetable(polygon_distribution_T, fullfile(outputDir, 'polygon_distribution.xls'))

        save(fullfile(outputDir,  'verticesInfo.mat'), 'cylindre2DImage', 'newVerticesNeighs2D', 'newVertices2D', 'centroids', 'validCellsFinal', 'borderCells', 'surfaceRatio', 'outputDir', 'nameOfSimulation', 'areaOfValidCells');
        save(fullfile(outputDir, 'allInfo.mat'), 'deployedImg', 'deployedImg3x', 'wholeImage', 'polygon_distribution', 'newNeighbours2D', 'newNeighbours2D_Checked');
    end
end

