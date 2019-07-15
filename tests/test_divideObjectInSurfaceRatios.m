function test_divideObjectInSurfaceRatios(testCase)
%TEST_1 Summary of this function goes here
%   Detailed explanation goes here
warning('off','all')

infoPerSurfaceRatio = divideObjectInSurfaceRatios('D:\Pablo\LimeSeg_Pipeline\data\Salivary gland_ExtractedVertices_Correct\Wildtype\2017-12-04\1a\Results\', 1);
actSolution = cell2table(infoPerSurfaceRatio, 'VariableNames', {'Image3DWithVolumen', 'SR3D', 'Layer3D', 'ApicalBasalCellFeatures3D', 'BasalApicalCellFeatures3D'});
load('D:\Pablo\LimeSeg_Pipeline\data\Salivary gland_ExtractedVertices_Correct\Wildtype\2017-12-04\1a\Results\glandDividedInSurfaceRatios_AllUnrollFeatures.mat', 'infoPerSurfaceRatio')
expSolution = infoPerSurfaceRatio(:, 1:5);
verifyEqual(testCase,actSolution,expSolution);

warning('on','all')
end

