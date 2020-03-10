function [] = showSelectedCell()
%SHOWSELECTEDCELL Summary of this function goes here
%   Detailed explanation goes here
selectCellId = getappdata(0, 'cellId');
labelledImage = getappdata(0, 'labelledImageTemp_Resized');
selectedZ = getappdata(0, 'selectedZ');
lumenImage = getappdata(0, 'lumenImage_Resized');

% perimImg = bwperim(labelledImage(:, :,  selectedZ) == selectCellId)';
%imshow(perimImg);
imageSequence = getappdata(0, 'imageSequence');

imgToShow = imageSequence(:, :, selectedZ)';

% imgToShow(perimImg == 1) = 65536;
cla('reset') 

imshow(uint16(imgToShow));
hold on;
    if selectCellId > 0
        [xIndices, yIndices] = find(labelledImage(:, :,  selectedZ) == selectCellId);
        if isempty(xIndices) == 0
            s2 = scatter(xIndices, yIndices, 'blue','filled','SizeData',10);
            hold off
            alpha(s2,.4)
        end
        [xIndices, yIndices] = find(lumenImage(:, :,  selectedZ) == 1);
    else
        [xIndices, yIndices] = find(lumenImage(:, :,  selectedZ) == 1);
    end


%% Showing lumen
if isempty(xIndices) == 0 && getappdata(0, 'hideLumen') == 0
    hold on
    s = scatter(xIndices, yIndices, 'red', 'filled','SizeData',10);
    hold off
    alpha(s,.5)
end

datacursormode 'off';
dcm_obj = datacursormode();
set(dcm_obj,'UpdateFcn',@pickCell);

end

