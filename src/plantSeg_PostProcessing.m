function plantSeg_PostProcessing(selpath)
%PLANTSEG_POSTPROCESSING Summary of this function goes here
%   Detailed explanation goes here

    selpath = 'D:\Pablo\LimeSeg_Pipeline\data\Cysts\cystPlantSeg\cyst10_predictions_gasp_average.tiff';
    
    tiff_info = imfinfo(selpath); % return tiff structure, one element per image
    tiff_stack = imread(selpath, 1) ; % read in first image
    %concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_info, 1)
        temp_tiff = imread(selpath, ii);
        tiff_stack = cat(3 , tiff_stack, temp_tiff);
    end
    image_labelled = double(tiff_stack)-1;
    
    lumenImage = image_labelled == 12; 
    image_labelled(image_labelled == 12) = 0;
    [labelledImage, basalLayer, apicalLayer] = postprocessGland(image_labelled, double(tiff_stack) == 1, lumenImage, '', colorcube(100), 0);
    
    calculateMissingCells(labelledImage, lumenImage, apicalLayer, basalLayer, colorcube(100), []);
end

