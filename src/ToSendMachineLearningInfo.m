%To Adrian and Maria Jose's Machine Learning project 
clear all

files = dir('**/data/Salivary gland/**/Results/3d_layers_info.mat');
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0 && (contains(lower(x), 'wildtype')), {files.folder}); %% contains(lower(x), 'echnoid')
files = files(nonDiscardedFiles);

for numFile = 1:length(files)
    inputDir = files(numFile).folder;
    inputDirSplitted = strsplit(inputDir, '\');
    if exist(fullfile(inputDirSplitted{1:end-1}, 'RawSequence')) == 0
        continue
    end
    inputDir = fullfile(inputDirSplitted{1:end-1});

    load(fullfile(inputDir, 'Results\3d_layers_info.mat'), 'labelledImage');
    labelledImage_NoTips = labelledImage(6:end-5, 6:end-5, 6:end-5);
    usedZScale = 1024/size(labelledImage_NoTips, 1)
    labelledImage3D = imresize3(labelledImage_NoTips, usedZScale, 'nearest');
    labelledImageAsConfocal = imresize3(labelledImage_NoTips, [1024 1024 size(labelledImage_NoTips, 3)], 'nearest');

    imageSequenceFiles = dir(fullfile(inputDir, 'RawSequence/*.tif'));
    NoValidFiles = startsWith({imageSequenceFiles.name},'._','IgnoreCase',true);
    imageSequenceFiles = imageSequenceFiles(~NoValidFiles);
    imageSequence = [];
    imageSequence_Small = [];

    for numImg = 1:size(imageSequenceFiles, 1)
        actualFile = imageSequenceFiles(numImg);
        imageSequence(:, :, numImg) =  imread(fullfile(actualFile.folder, actualFile.name))';
        imageSequence_Small(:, :, numImg) =  imresize(imread(fullfile(actualFile.folder, actualFile.name))', [size(labelledImage_NoTips, 1) size(labelledImage_NoTips, 2)], 'nearest');
    end

    rois = ReadImageJROI(fullfile(inputDir, 'Cells\RoiSetCells.zip'));
    %vnRectBounds = ['nTop', 'nLeft', 'nBottom', 'nRight']
    for numRoi = 1:length(rois)
        centroidOfRois(numRoi, 1:3) = [round(mean(rois{numRoi}.vnRectBounds([1 3]))) round(mean(rois{numRoi}.vnRectBounds([2 4]))) rois{numRoi}.nPosition];
        centroidOfRois(numRoi, 4) = abs(centroidOfRois(numRoi, 1) - rois{numRoi}.vnRectBounds(1));
    %     p = figure, imshow(double(imageSequence(:, :, centroidOfRois(numRoi, 3)))/65535)
    %     hold on;
    %     plot(centroidOfRois(numRoi, 1), centroidOfRois(numRoi, 2), 'rx');
    %     imline(gca, [centroidOfRois(numRoi, 1:2)+centroidOfRois(numRoi, 4); centroidOfRois(numRoi, 1:2)])
    end
    
    centroidOfRois = array2table(centroidOfRois, 'VariableNames', {'CoordX', 'CoordY', 'CoordZ', 'RadiusSphere'});
    imageSequenceInterpolated = round(imresize3(imageSequence, [1024 1024 size(labelledImage3D, 3)], 'linear'));
    
    save(fullfile(inputDir, strcat(strjoin(inputDirSplitted(end-3:end-1), '_'),'.mat')), 'labelledImage3D', 'imageSequence', 'imageSequenceInterpolated', 'centroidOfRois', 'usedZScale', '-v7.3');

%     zFrame = 35;
%     %% Resized img
%     figure, imshow(double(imageSequence(:, :, zFrame))/4095)
%     hold on;
%     [xIndices, yIndices] = find(labelledImageAsConfocal(:, :,  zFrame)>0);
%     if isempty(xIndices) == 0
%         s2 = scatter(yIndices, xIndices, 'blue','filled','SizeData',10);
%         hold off
%         alpha(s2,.1)
%     end

    % %% Small img
    % figure, imshow(double(imageSequence_Small(:, :, zFrame))/65535)
    % hold on;
    % [xIndices, yIndices] = find(labelledImage_NoTips(:, :,  zFrame)>0);
    % if isempty(xIndices) == 0
    %     s2 = scatter(yIndices, xIndices, 'blue','filled','SizeData',10);
    %     hold off
    %     alpha(s2,.1)
    % end
end
