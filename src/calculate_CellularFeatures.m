function [CellularFeaturesWithNoValidCells, meanSurfaceRatio] = calculate_CellularFeatures(neighbours_data,apical3dInfo,basal3dInfo,apicalLayer,basalLayer,labelledImage,noValidCells,validCells,polygon_distribution,outputDir)
%CALCULATE_CELLULARFEATURES Summary of this function goes here
%   Detailed explanation goes here
%%  Calculate number of neighbours of each cell
number_neighbours=table(cellfun(@length,(apical3dInfo.neighbourhood)),cellfun(@length,(basal3dInfo.neighbourhood)));
total_neighbours3D=calculateNeighbours3D(labelledImage);
total_neighbours3DRecount=cellfun(@(x) length(x), total_neighbours3D.neighbourhood, 'UniformOutput',false);
apicobasal_neighbours=cellfun(@(x,y)(unique(vertcat(x,y))), apical3dInfo.neighbourhood, basal3dInfo.neighbourhood, 'UniformOutput',false);
apicobasal_neighboursRecount=cellfun(@(x) length(x),apicobasal_neighbours,'UniformOutput',false);

%%  Calculate area cells
apical_area_cells=cell2mat(struct2cell(regionprops(apicalLayer,'Area'))).';
basal_area_cells=cell2mat(struct2cell(regionprops(basalLayer,'Area'))).';
surfaceRatio = basal_area_cells ./ apical_area_cells;
%meanSurfaceRatio = mean(surfaceRatioValidCells);
meanSurfaceRatio = sum(basal_area_cells(validCells)) / sum(apical_area_cells(validCells));

%%  Calculate volume cells
volume_cells=table2array(regionprops3(labelledImage,'Volume'));

%%  Determine if a cell is a scutoid or not
scutoids_cells=cellfun(@(x,y) double(~isequal(x,y)), neighbours_data.Apical,neighbours_data.Basal);

%%  Export to a excel file
ID_cells=(1:length(basal3dInfo.neighbourhood)).';

 if isequal(total_neighbours3D.neighbourhood,apicobasal_neighbours)==0
        
        pos=cellfun(@isequal, total_neighbours3D.neighbourhood,apicobasal_neighbours);
       
        ids=ID_cells(pos==0);
        ids(ismember(ids,noValidCells))=[];
        
        
        IDsStrings=string(num2str(ids));
        IDsStrings=strjoin(IDsStrings,', ');
        
        msg1="Cells with IDs ";
        msg2=strcat(msg1,IDsStrings);
     
        msg3="  could be wrong due to Total_neighbours is different from Apicobasal_neighours";
        msg=strcat(msg2,msg3);
    
        warning(msg);
 end

CellularFeatures=table(ID_cells,number_neighbours.Var1,number_neighbours.Var2,total_neighbours3DRecount,apicobasal_neighboursRecount,scutoids_cells,apical_area_cells,basal_area_cells, surfaceRatio, volume_cells);
CellularFeatures.Properties.VariableNames = {'ID_Cell','Apical_sides','Basal_sides','Total_neighbours','Apicobasal_neighbours','Scutoids','Apical_area','Basal_area', 'Surface_Ratio','Volume'};
CellularFeaturesWithNoValidCells = CellularFeatures;
CellularFeatures(noValidCells,:)=[];


if isempty(outputDir) == 0
    writetable(CellularFeatures,fullfile(outputDir,'Results', 'cellular_features_LimeSeg3DSegmentation.xls'), 'Range','B2');

    %% Poligon distribution 
    polygon_distribution_3D=calculate_polygon_distribution(cellfun(@length, total_neighbours3D.neighbourhood), validCells);
    writetable(table('','VariableNames',{'Apical'}),fullfile(outputDir,'Results', 'cellular_features_LimeSeg3DSegmentation.xls'), 'Sheet', 2, 'Range', 'B2')
    writetable(table(polygon_distribution.Apical),fullfile(outputDir,'Results', 'cellular_features_LimeSeg3DSegmentation.xls'), 'Sheet', 2, 'Range', 'B3', 'WriteVariableNames',false);
    writetable(table('','VariableNames',{'Basal'}),fullfile(outputDir,'Results', 'cellular_features_LimeSeg3DSegmentation.xls'), 'Sheet', 2, 'Range', 'B6')
    writetable(table(polygon_distribution.Basal),fullfile(outputDir,'Results', 'cellular_features_LimeSeg3DSegmentation.xls'), 'Sheet', 2, 'Range', 'B7', 'WriteVariableNames',false);
    writetable(table('','VariableNames',{'Accumulate'}),fullfile(outputDir,'Results', 'cellular_features_LimeSeg3DSegmentation.xls'), 'Sheet', 2, 'Range', 'B10')
    writetable(table(polygon_distribution_3D),fullfile(outputDir,'Results', 'cellular_features_LimeSeg3DSegmentation.xls'), 'Sheet', 2, 'Range', 'B11', 'WriteVariableNames',false);
end
