path2load = dir('D:\Pablo\LimeSeg_Pipeline\data\Salivary gland\Wildtype\**\Wildtype*.mat');
path2save = 'C:\Users\PhD-Student\Desktop\imagesStarDist\';

for numFile = 1:size(path2load,1)
   
    load([path2load(numFile).folder '\' path2load(numFile).name],'imageSequence','imageSequenceInterpolated','labelledImage3D')
    
    [z1]=size(imageSequence,3);
    [z2]=size(imageSequenceInterpolated,3);
    [z3]=size(labelledImage3D,3);

    stepZ = round(z3/z1);
    
    imageSequenceInterpolatedStepZ = uint16(imageSequenceInterpolated(:,:,1:stepZ:end));
    [H,W,c]=size(imageSequenceInterpolatedStepZ);
    labelledImage3DStepZ = imresize3(uint16(labelledImage3D(:,:,1:stepZ:end)),[H,W,c],'nearest');
  
    %write a Tiff file, appending each image as a new page
    for ii = 1 : size(imageSequenceInterpolatedStepZ, 3)
        imwrite(imageSequenceInterpolatedStepZ(:,:,ii) , [path2save 'raw\rawImage' num2str(numFile) '.tif'] , 'WriteMode' , 'append') ;
        imwrite(labelledImage3DStepZ(:,:,ii) , [path2save 'mask\labelImage' num2str(numFile) '.tif'] , 'WriteMode' , 'append') ;
    end
end

