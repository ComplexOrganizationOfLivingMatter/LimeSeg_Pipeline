function [] = showSelectedCell()
%SHOWSELECTEDCELL Summary of this function goes here
%   Detailed explanation goes here
selectCellId = getappdata(0, 'cellId');
labelledImage = getappdata(0, 'labelledImageTemp_Resized');
selectedZ = getappdata(0, 'selectedZ');
lumenImage = getappdata(0, 'lumenImage_Resized');
showAllCells = getappdata(0, 'showAllCells');
imageSequence = getappdata(0, 'imageSequence');
colours = getappdata(0, 'colours');

imgToShow = mat2gray(imageSequence(:, :, selectedZ)');

cla

if showAllCells==1
    %% Showing all cells
    labImageZ = labelledImage(:, :,  selectedZ)';
    centLab = cat(1,regionprops(labImageZ,'Centroid'));
    centroids = vertcat(centLab.Centroid);
    labelsZ = unique(labImageZ);
    
    B = labeloverlay(imgToShow,labelledImage(:, :,  selectedZ)', 'Colormap', colours);
    imshow(B);
    hold on;
    textscatter(centroids(labelsZ(2:end),1),centroids(labelsZ(2:end),2),cellfun(@num2str,num2cell(labelsZ(2:end)),'UniformOutput',false),'TextDensityPercentage',100,'ColorData',ones(length(labelsZ(2:end)),3));
else
    imshow(imgToShow);
    if selectCellId > 0
        [xIndices, yIndices] = find(labelledImage(:, :,  selectedZ) == selectCellId);
        if isempty(xIndices) == 0
            hold on
            s2 = scatter(xIndices, yIndices, 'blue','filled','SizeData',10);
            hold off
            alpha(s2,.4)
        end
    end
    hold off
end


%% Showing lumen
[xIndices, yIndices] = find(lumenImage(:, :,  selectedZ) == 1);
if isempty(xIndices) == 0 && getappdata(0, 'hideLumen') == 0
    hold on
    s = scatter(xIndices, yIndices, 'red', 'filled','SizeData',10);
    hold off
    alpha(s,.5)
end


end

