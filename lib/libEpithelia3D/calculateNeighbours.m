function [neighs_real,sides_cells]=calculateNeighbours(L_img,varargin)

%% Generate neighbours
    ratio=4;
    se = strel('disk',ratio);
    neighs_real={};
    cells=sort(unique(L_img));
    cells=cells(cells>0);                  %% Deleting cell 0 from range

    if ~isempty(varargin)
        se = strel('disk',varargin{1});
    end

    for cel = cells'
        BW = bwperim(L_img==cel);
        BW_dilate=imdilate(BW,se);
        neighs=unique(L_img(BW_dilate==1));
        neighs_real{cel}=neighs((neighs ~= 0 & neighs ~= cel));
        sides_cells(cel)=length(neighs_real{1,cel});

    end


end