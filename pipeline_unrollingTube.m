%% Unroll tube

addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))
%addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files = dir('**/Salivary gland/**/Results/3d_layers_info.mat');

% for numFile = 1:length(files)
%     if contains(lower(files(numFile).folder), 'discarded') == 0
%         load(fullfile(files(numFile).folder, '3d_layers_info.mat'));
% 
%         figure; paint3D(basalLayer, [], colours)
%         xlabel('x')
%         ylabel('y')
%         zlabel('z')
%         axisToChange = input('axis:');
%         if isequal(axisToChange, '') == 0
%             
%             size(basalLayer)
%             valuesToChange = input('values to remove:');
%             if isequal(axisToChange, 'x')
%                 basalLayer(valuesToChange, :, :) = 0;
%                 apicalLayer(valuesToChange, :, :) = 0;
%                 labelledImage(valuesToChange, :, :) = 0;
%                 lumenImage(valuesToChange, :, :) = 0;
%             elseif isequal(axisToChange, 'y')
%                 basalLayer(:, valuesToChange, :) = 0;
%                 apicalLayer(:, valuesToChange, :) = 0;
%                 labelledImage(:, valuesToChange, :) = 0;
%                 lumenImage(:, valuesToChange, :) = 0;
%             end
%             
%             figure; paint3D(basalLayer, [], colours)
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
% 
%             resultOK = input('everithing ok? ');
%             if isequal(resultOK, 'OK')
%                 save(fullfile(files(numFile).folder, '3d_layers_info.mat'), 'labelledImage', 'basalLayer', 'apicalLayer', 'apical3dInfo', 'basal3dInfo', 'colours', 'lumenImage','glandOrientation', '-v7.3')
%             end
%         end
%     end
%     close all
% end

nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0 && contains(lower(x), 'wildtype'), {files.folder});
files = files(nonDiscardedFiles);

parfor numFile = 1:length(files)
    files(numFile).folder
    selpath = files(numFile).folder;
    
    unrollTube_parallel(selpath);
end