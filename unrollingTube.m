%% Unroll tube

addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))
%addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files = dir('**/Salivary gland/**/Results/3d_layers_info.mat');

for numFile = 1:length(files)
    if contains(lower(files(numFile).folder), 'discarded') == 0
        files(numFile).folder
        selpath = files(numFile).folder;
        
        unrollTube_parallel(selpath);
    end
end