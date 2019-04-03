function [polygon_distribution, neighbours_data] = limeSeg_PostProcessing(outputDir)
%PIPELINE Summary of this function goes here
%   Detailed explanation goes here
    mkdir(fullfile(outputDir, 'Cells', 'OutputLimeSeg'));
    mkdir(fullfile(outputDir, 'ImageSequence'));
    mkdir(fullfile(outputDir, 'Lumen', 'SegmentedLumen'));
    mkdir(fullfile(outputDir, 'Results'));
    mkdir(fullfile(outputDir, 'Apical_Labelled'));

    zScale = inputdlg('Insert z-scale of Gland');
    zScale = zScale{1};
    
    save(fullfile(outputDir, 'zScaleOfGland.mat'), 'zScale');
    
    resizeImg = 1/zScale;

    tipValue = 4;

    imageSequenceFiles = dir(fullfile(outputDir, 'ImageSequence/*.tif'));
    NoValidFiles = startsWith({imageSequenceFiles.name},'._','IgnoreCase',true);
    imageSequenceFiles=imageSequenceFiles(~NoValidFiles);
    demoFile =  imageSequenceFiles(3);
    demoImg = imread(fullfile(demoFile.folder, demoFile.name));

    imgSize = round(size(demoImg)*resizeImg);

    if exist(fullfile(outputDir, 'Results', '3d_layers_info.mat'), 'file')
        load(fullfile(outputDir, 'Results', '3d_layers_info.mat'))
    else
        colours = [];
        [labelledImage, outsideGland] = processCells(fullfile(outputDir, 'Cells', filesep), resizeImg, imgSize, tipValue);

        [labelledImage, lumenImage] = processLumen(fullfile(outputDir, 'Lumen', filesep), labelledImage, resizeImg, tipValue);
        labelledImage = fill0sWithCells(labelledImage, imclose(labelledImage, strel('sphere', 3)) == 0);
            
        %% Put both lumen and labelled image at a 90 degrees

        orientationGland = regionprops3(lumenImage>0, 'Orientation');
        glandOrientation = -orientationGland.Orientation(1);
        %labelledImage = imrotate(labelledImage, glandOrientation);
        %lumenImage = imrotate(lumenImage, glandOrientation);

        %% Get basal layer by dilating the empty space
        [basalLayer] = getBasalFrom3DImage(labelledImage, lumenImage, tipValue);

        %% Get apical layer by dilating the lumen
        [apicalLayer] = getApicalFrom3DImage(lumenImage, labelledImage);
        exportAsImageSequence(apicalLayer, fullfile(outputDir, 'Apical_Labelled'), colours, tipValue);

        %% Export image sequence
        [colours] = exportAsImageSequence(labelledImage, fullfile(outputDir, 'Cells', 'labelledSequence', filesep), colours, tipValue);
    end
    [outsideGland] = getOutsideGland(labelledImage);

    setappdata(0,'outputDir', outputDir);
    setappdata(0,'labelledImage',labelledImage);
    setappdata(0,'lumenImage', lumenImage);
    setappdata(0,'resizeImg',resizeImg);
    setappdata(0,'tipValue', tipValue);
    setappdata(0, 'glandOrientation', glandOrientation);

    if exist(fullfile(outputDir, 'Results', 'valid_cells.mat'), 'file')
        load(fullfile(outputDir, 'Results', 'valid_cells.mat'))
    else
        [noValidCells] = insertNoValidCells();
        validCells = setdiff(1:max(labelledImage(:)), noValidCells);
        save(fullfile(outputDir, 'Results', 'valid_cells.mat'), 'noValidCells', 'validCells')
    end
    [answer, apical3dInfo, notFoundCellsApical, basal3dInfo, notFoundCellsBasal] = calculateMissingCells(labelledImage, lumenImage, apicalLayer, basalLayer, colours, noValidCells);

    %% Insert no valid cells
    while isequal(answer, 'Yes')
        h = window();
        waitfor(h);

        savingResults = saveResults();

        if isequal(savingResults, 'Yes')
            labelledImage = getappdata(0, 'labelledImageTemp');
            close all
            [labelledImage] = fillEmptySpacesByWatershed3D(labelledImage, outsideGland | lumenImage, 1);
            outsideGland_NotLumen = ~outsideGland | imdilate(lumenImage, strel('sphere', 2));
            
            labelledImage = fill0sWithCells(labelledImage, outsideGland | lumenImage);
            labelledImage = fill0sWithCells(labelledImage, outsideGland_NotLumen == 0);
            labelledImage(lumenImage) = 0;
            exportAsImageSequence(labelledImage, fullfile(outputDir, 'Cells', 'labelledSequence', filesep), colours, tipValue);

            %% Calculate neighbours and plot missing cells
            [basalLayer] = getBasalFrom3DImage(labelledImage, lumenImage, tipValue);
            [apicalLayer] = getApicalFrom3DImage(lumenImage, labelledImage);
            exportAsImageSequence(apicalLayer, fullfile(outputDir, 'Apical_Labelled'), colours, tipValue);
            [answer, apical3dInfo, notFoundCellsApical, basal3dInfo, notFoundCellsBasal] = calculateMissingCells(labelledImage, lumenImage, apicalLayer, basalLayer, colours, noValidCells);
        else
            [answer] = isEverythingCorrect();
        end
        setappdata(0,'labelledImage',labelledImage);
    end

    %% Save apical and basal 3d information
    save(fullfile(outputDir, 'Results', '3d_layers_info.mat'), 'labelledImage', 'basalLayer', 'apicalLayer', 'apical3dInfo', 'basal3dInfo', 'colours', 'lumenImage','glandOrientation', '-v7.3')

    %% Calculate poligon distribution
    [polygon_distribution_Apical] = calculate_polygon_distribution(cellfun(@length, apical3dInfo), validCells);
    [polygon_distribution_Basal] = calculate_polygon_distribution(cellfun(@length, basal3dInfo), validCells);
    neighbours_data = table(apical3dInfo, basal3dInfo);
    polygon_distribution = table(polygon_distribution_Apical, polygon_distribution_Basal);
    neighbours_data.Properties.VariableNames = {'Apical','Basal'};
    polygon_distribution.Properties.VariableNames = {'Apical','Basal'};

    %% Export to excel cellular features
    cellularFeatures = calculate_CellularFeatures(neighbours_data,apical3dInfo,basal3dInfo,apicalLayer,basalLayer,labelledImage,noValidCells,validCells,polygon_distribution,outputDir);
    
%   save(fullfile(outputDir, 'Results', 'cellularFeaturesExcel.mat'), cellularFeatures); 
end

