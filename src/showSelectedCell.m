function [] = showSelectedCell()
%SHOWSELECTEDCELL Summary of this function goes here
%   Detailed explanation goes here
selectCellId = getappdata(0, 'cellId');
labelledImage = getappdata(0, 'labelledImageTemp_Resized');
selectedZ = getappdata(0, 'selectedZ');
lumenImage = getappdata(0, 'lumenImage_Resized');
showAllCells = getappdata(0, 'showAllCells');
cmap = getappdata(0, 'cmap');
% perimImg = bwperim(labelledImage(:, :,  selectedZ) == selectCellId)';
%imshow(perimImg);

imageSequence = getappdata(0, 'imageSequence');
imgToShow = imageSequence(:, :, selectedZ)';

% imgToShow(perimImg == 1) = 65536;
cla
%clf

if showAllCells==1
    %% Showing all cells
    labImageZ = labelledImage(:, :,  selectedZ)';
    centLab = cat(1,regionprops(labImageZ,'Centroid'));
    centroids = vertcat(centLab.Centroid);
    labelsZ = unique(labImageZ);
    
    cmap(1,:)=[0 0 0];
    image(labelledImage(:, :,  selectedZ)');
%     hImage1 = imshow(uint16(imgToShow)); 
    colormap(cmap)

    %imshow(labelledImage(:, :,  selectedZ)',cmap);
    hold on
%     hImage2 = imshow(labelledImage(:, :,  selectedZ)',cmap); 
%     set(hImage2, 'AlphaData', 0.35);
    im = image(imgToShow);
    im.AlphaData = 0.4;
    

    textscatter(centroids(labelsZ(2:end),1),centroids(labelsZ(2:end),2),cellfun(@num2str,num2cell(labelsZ(2:end)),'UniformOutput',false),'TextDensityPercentage',100,'ColorData',ones(length(labelsZ(2:end)),3));
    hold off

else
    image(imgToShow); 
    colormap('gray');
    hold on
    if selectCellId > 0
        [xIndices, yIndices] = find(labelledImage(:, :,  selectedZ) == selectCellId);
        if isempty(xIndices) == 0
            hold on
            s2 = scatter(xIndices, yIndices, 'blue','filled','SizeData',10);
            hold off
            alpha(s2,.4)
        end
    end
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

