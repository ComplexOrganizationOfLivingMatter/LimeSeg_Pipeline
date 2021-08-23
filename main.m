%% First pipeline to modify mistakes on S. Glands
addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))

close all

selpath = uigetdir('data');
if isempty(selpath) == 0
    limeSeg_PostProcessing(selpath);
end

