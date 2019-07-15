function tests = executesTests
    addpath(genpath('tests'));
    localfunctions = matlab.unittest.TestCase.forInteractiveUse;
    tests = test_1(localfunctions);
end