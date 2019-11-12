function exportLumen(lumenImage, outputDir, tipValue)
%EXPORTASIMAGESEQUENCE Summary of this function goes here
%   Detailed explanation goes here

    mkdir(fullfile(outputDir, 'Lumen', 'inferLumen'));        
    
    h = figure('Visible', 'off');
    for numZ = 1+tipValue+1:(size(lumenImage, 3)-(tipValue+1))
        imshow(imcomplement((lumenImage(:, :, numZ))));
        set(h, 'units','normalized','outerposition',[0 0 1 1]);
        ax = get(h, 'Children');
        set(ax,'Units','normalized')
        set(ax,'Position',[0 0 1 1])
        h.InvertHardcopy = 'off';
        saveas(h,fullfile(outputDir,'Lumen\inferLumen\', strcat('lumenImage_', num2str(numZ-(tipValue+1)), '.tif')))
        %imwrite(lumenImage(:, :, numZ), , );
    end
    
end