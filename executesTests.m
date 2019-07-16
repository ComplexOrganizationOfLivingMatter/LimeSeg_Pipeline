function tests = executesTests
    addpath(genpath('tests'));
    mode = matlab.unittest.TestCase.forInteractiveUse;
    
    perfectSample = 'D:\Pablo\LimeSeg_Pipeline\data\Salivary gland_ExtractedVertices_Correct\Wildtype\2017-12-04\1a\Results\';
    
    %Testing function 'divideObjectInSurfaceRatios'
    %test_divideObjectInSurfaceRatios(mode);
    
    %% Testing function 'mappCylindricalCoordinatesInto2D' from 'unrollTube'
    % Apical
    test_mapped3Dinto2D_apical(mode, perfectSample, '1')
    test_mapped3Dinto2D_apical(mode, perfectSample, 'basal')
    
end