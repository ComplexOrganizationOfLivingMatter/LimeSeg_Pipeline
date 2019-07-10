function img3d = resizeEachCell(img3d_original, img3d, img3dComplete, imgSize)
%RESIZEEACHCELL Summary of this function goes here
%   Detailed explanation goes here

    newPixelsPerCell = cell(max(img3d_original(:)), 1);
    for numCell = 1:max(img3d_original(:))
        actualCellSmall = double(img3d_original == numCell);
        
        actualCellBig = imresize3(actualCellSmall, imgSize, 'Method', 'nearest');
        
        newPixelsPerCell{numCell} = find(img3d == 0 & actualCellBig>0);
    end
    
    %% See if any new pixel of a cell intersect with any other cell
    for numCell = 1:max(img3d_original(:))
        newPixelsInteresect{numCell} = find(cellfun(@(x) isempty(intersect(newPixelsPerCell{numCell}, x))==0, newPixelsPerCell));
    end
    newPixelsInteresect
end

