%% First pipeline to modify mistakes on S. Glands
addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))

close all

[fileName, selpath] = uigetfile('*.*');
if isempty(selpath) == 0
    limeSeg_PostProcessing(selpath, fileName);
end

