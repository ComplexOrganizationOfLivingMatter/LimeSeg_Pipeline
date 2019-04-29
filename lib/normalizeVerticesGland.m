function [samiraTable, glandMinValue] = normalizeVerticesGland(samiraTable, cutValue)
%NORMALIZEVERTICESGLAND Summary of this function goes here
%   Detailed explanation goes here

for numRow=1:size(samiraTable,1)
    actualVertices = samiraTable{numRow,5};
    actualVertices(1:2:end)= actualVertices(1:2:end) - cutValue;
    samiraTable{numRow,5}=actualVertices;
    minValues(numRow)= min(actualVertices);
    if cutValue < 0 & minValues(numRow) < 0 
        disp('error');
    end
end 

glandMinValue=min(minValues);

end

