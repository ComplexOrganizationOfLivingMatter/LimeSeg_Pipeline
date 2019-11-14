function exportLumen(lumenImage, outputDir, tipValue)
%EXPORTASIMAGESEQUENCE Summary of this function goes here
%   Detailed explanation goes here

    mkdir(fullfile(outputDir, 'Lumen', 'inferLumen'));        
  
    for numZ = 1+tipValue+1:(size(lumenImage, 3)-(tipValue+1))
        actualImg = imcomplement(imresize(lumenImage(:, :, numZ)',[1024,1024], 'nearest'));
        imwrite(double(actualImg), fullfile(outputDir,'Lumen\inferLumen\', strcat('lumenImage_', num2str(numZ-(tipValue+1)), '.tif')))
    end
    
end