function [rotatedImg3d, rotations] = rotateImg3(img3d, apicalRotationsOriginal)

    rotations = [0 0 0];
    closedImg3D = imdilate(img3d>0,strel('sphere',5));
    orientationImg3D = imerode(bwmorph3(closedImg3D,'fill'),strel('sphere',3));
    clearvars closedImg3D
    
    %rotation A %% You may keep this untouched
    orientationObj = regionprops3(orientationImg3D, 'Orientation');
    if exist('apicalRotationsOriginal', 'var') == 0
        img3d_rot1 = imrotate(img3d, - orientationObj.Orientation(1));
        orientationImg3D = imrotate(img3d, - orientationObj.Orientation(1));
    else
        img3d_rot1 = imrotate(img3d, - apicalRotationsOriginal(1));
        orientationImg3D = imrotate(img3d, - apicalRotationsOriginal(1));
    end
    rotations(1) = orientationObj.Orientation(1);
    
    clearvars img3d
    
    %rotation B
    xzyImg = permute(img3d_rot1,[1 3 2]);
    orientationImg3D = permute(orientationImg3D,[1 3 2]);
    orientationObj = regionprops3(orientationImg3D, 'Orientation');
    if exist('apicalRotationsOriginal', 'var') == 0
        xzyImg3d_rot = imrotate(xzyImg, - orientationObj.Orientation(1));
    else
        xzyImg3d_rot = imrotate(xzyImg, - apicalRotationsOriginal(2));
    end
    rotations(2) = orientationObj.Orientation(1);
    clearvars xzyImg img3d_rot1
    
    %rotation C
    yzxImg3d_rot = permute(xzyImg3d_rot,[3 2 1]);
    orientationImg3D = permute(orientationImg3D,[3 2 1]);

    orientationObj = regionprops3(orientationImg3D, 'Orientation');
    if exist('apicalRotationsOriginal', 'var') == 0
        img3d_rotFinal = imrotate(yzxImg3d_rot, - orientationObj.Orientation(1));
    else
        img3d_rotFinal = imrotate(yzxImg3d_rot, - apicalRotationsOriginal(3)); 
    end
    rotations(3) = orientationObj.Orientation(1);
    clearvars xzyImg3d_rot orientationImg3D
    
    %come back to original axes
    img3d_rot = permute(img3d_rotFinal,[3 1 2]);
    clearvars img3d_rotFinal
    
    [xCol, yRow, z] = ind2sub(size(img3d_rot),find(img3d_rot>0));
    if min(xCol) <= 0 || min(yRow) <= 0 || min(z) <= 0

        if min(xCol) <= 0 
            xCol = xCol + abs(min(xCol)) + 1;
        end

        if min(yRow) <= 0
            yRow = yRow + abs(min(yRow)) + 1;
        end

        if min(z) <= 0
            z = z + abs(min(z)) +1;
        end
        rotatedImg3d = zeros(max(xCol),max(yRow),max(z));
        rotatedImg3d( sub2ind(size(rotatedImg3d),xCol,yRow,z) ) = img3d_rot( img3d_rot > 0 );
    else
        rotatedImg3d = img3d_rot;
    end
    
end