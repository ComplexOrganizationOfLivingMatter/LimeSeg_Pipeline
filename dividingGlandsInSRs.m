addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))
%addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files = dir('**/Salivary gland/**/Results/3d_layers_info.mat');

parfor numFile = 1:length(files)
    if contains(lower(files(numFile).folder), 'discarded') == 0
        
        selpath = files(numFile).folder;
        divideObjectInSurfaceRatios(selpath);
    end
end

