function [samiraTable] = normalizeVerticesGland(samiraTable,deployedImg3x)
%NORMALIZEVERTICESGLAND Summary of this function goes here
%   Detailed explanation goes here
for numRow=1:size(samiraTable,1)
    actualVertices = samiraTable{numRow,5};
    actualVertices(1:2:end)= actualVertices(1:2:end)- round(size(deployedImg3x,1)/2);
    samiraTable{numRow,5}=actualVertices;
    minValues(numRow)= min(actualVertices);
end 
glandMinValue=min(minValues);

end

