clear all
modelToLoad = {'E-cadh Inhibited/e-cadhi type I (flatten severe)','E-cadh Inhibited/e-cadhi type II (flatten intermediate)','Echinoid Inhibited'};

% path2load = dir('D:\Pablo\LimeSeg_Pipeline\data\Salivary gland\Wildtype\**\Wildtype*.mat');
pathGlands = 'E:\Pedro\SalivaryGlands\';

tipValue = 5;
xySize = 1024;

for nTissue = 1:length(modelToLoad)
    path2load = dir(['../data/Salivary gland/' modelToLoad{nTissue} '/**/3d_layers_info.mat']);
    
    path2save = fullfile(pathGlands,modelToLoad{nTissue});
    mkdir(fullfile(path2save,'labelledImage'))
    mkdir(fullfile(path2save,'lumenImage'))
    
    if contains(modelToLoad{nTissue},'severe') || contains(modelToLoad{nTissue},'intermediate') || contains(modelToLoad{nTissue},'Echinoid')
        selectedIDs = ~cellfun(@(x) contains(lower(x),'discarded'),{path2load(:).folder});
        
        if contains(modelToLoad{nTissue},'intermediate')
            selectedIDs = cellfun(@(x) contains(lower(x),'oldmethod'),{path2load(:).folder});
        end
        path2load_1 = path2load(selectedIDs);
    
        
        for numFile = 1:size(path2load_1,1)

            load([path2load_1(numFile).folder '\' path2load_1(numFile).name],'labelledImage','lumenImage','glandOrientation')

            nameGland = strsplit(path2load_1(numFile).folder,'\');
            nameGland = [nameGland{end-2} '_' nameGland{end-1}];

            labelledImageWithoutTips = labelledImage(tipValue:end-tipValue,tipValue:end-tipValue,tipValue+1:end-tipValue);
            lumenWithoutTips = lumenImage(tipValue:end-tipValue,tipValue:end-tipValue,tipValue+1:end-tipValue);

            listImageSeq = dir([strrep(path2load_1(numFile).folder,'Results','ImageSequence/'),'*tif']);

            c = size(listImageSeq,1);

            labelledImageWithoutTips2 = fliplr(cat(3,labelledImageWithoutTips,zeros(size(labelledImageWithoutTips,1),size(labelledImageWithoutTips,2),c-size(labelledImageWithoutTips,3))));
            labelledImageWithoutTips2 = imrotate(labelledImageWithoutTips2,90);
            labelledImageResize = imresize3(uint16(labelledImageWithoutTips2),[xySize,xySize,c],'nearest');

            lumenWithoutTips2= fliplr(cat(3,lumenWithoutTips,zeros(size(lumenWithoutTips,1),size(lumenWithoutTips,2),c-size(lumenWithoutTips,3))));
            lumenWithoutTips2 = imrotate(lumenWithoutTips2,90);
            lumenImageResize = imresize3(uint16(lumenWithoutTips2),[xySize,xySize,c],'nearest');

            %write a Tiff file, appending each image as a new page
            for ii = 1 : size(labelledImageResize, 3)
                imwrite(labelledImageResize(:,:,ii) , fullfile(path2save, 'labelledImage',[nameGland '.tif']), 'WriteMode' , 'append') ;
                imwrite(lumenImageResize(:,:,ii) , fullfile(path2save, 'lumenImage',[nameGland '.tif']), 'WriteMode' , 'append') ;
            end
        end
    end

    if contains(modelToLoad{nTissue},'intermediate')
        selectedIDs = ~cellfun(@(x) contains(lower(x),'oldmethod'),{path2load(:).folder});
        path2load_2 = path2load(selectedIDs);
        
        for numFile = 1:size(path2load_2,1)
            load([path2load_2(numFile).folder '\' path2load_2(numFile).name],'labelledImage','lumenImage')
            nameGland = strsplit(path2load_2(numFile).folder,'\');
            nameGland = [nameGland{end-2} '_' nameGland{end-1}];
            %write a Tiff file, appending each image as a new page
            for ii = 1 : size(labelledImage, 3)
                imwrite(uint16(labelledImage(:,:,ii)) , fullfile(path2save, 'labelledImage',[nameGland '.tif']), 'WriteMode' , 'append') ;
                imwrite(uint16(lumenImage(:,:,ii)) , fullfile(path2save, 'lumenImage',[nameGland '.tif']), 'WriteMode' , 'append') ;
            end
            
        end
        
        
    end
    
    

end
