addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))

close all
clear all

files = dir('**/data/Salivary gland/**/global_3DFeatures.mat');

allFilesName = [];
allShapiroWilkHypothesis = table();
% allLeveneTests
globalMeanFeatures= {};
allTtestHypothesis = table();

for numFiles=1:length(files)
    
load(fullfile(files(numFiles).folder, 'global_3DFeatures.mat'), 'totalMeanFeatures'); 

totalMeanFeatures = removevars(totalMeanFeatures, {'Fun_PrincipalAxisLength'});
globalMeanFeatures{numFiles,1} = totalMeanFeatures;
featuresNames = fieldnames(totalMeanFeatures);

%% Shapiro-Wilk test to establish if the values follow a normal distribution.
featuresShapiroWilkHypothesis =table();
for nFeatures=1:size(totalMeanFeatures,2)
    [shapiroWilkHypothesis,shapiroWilkPvalue] = shapiroWilkTest(totalMeanFeatures{:,nFeatures});
    featuresShapiroWilkHypothesis= [featuresShapiroWilkHypothesis, mergevars(table(shapiroWilkHypothesis, shapiroWilkPvalue), [1 2], 'NewVariableName',featuresNames{nFeatures,1})];    
end

fileName = strsplit(files(numFiles).folder, {'/','\'});
fileName = convertCharsToStrings(fileName{end});
allFilesName = [allFilesName ; fileName];

allShapiroWilkHypothesis = [allShapiroWilkHypothesis ; featuresShapiroWilkHypothesis];

end

%% Student-t test to determine if the means of two populations are equal. 
% Only if there are not the same number of samples. Before we check if the
% samples share the same variance. If not, we apply the Welch correction to Student T
% test.

featuresVarianceHypothesis =table();

for nFeatures=1:size(totalMeanFeatures,2)

    if allShapiroWilkHypothesis{1,nFeatures}(1) == 0 && allShapiroWilkHypothesis{2,nFeatures}(1) == 0
        [varianceHypothesis,variancePvalue] = vartest2(globalMeanFeatures{1,1}{:,nFeatures},globalMeanFeatures{2,1}{:,nFeatures});
        featuresVarianceHypothesis= [featuresVarianceHypothesis, mergevars(table(varianceHypothesis, variancePvalue), [1 2], 'NewVariableName',featuresNames{nFeatures,1})];
        
        if varianceHypothesis == 0
        [tTestHypothesis,pValueTtest] = ttest2(globalMeanFeatures{1,1}{:,nFeatures},globalMeanFeatures{2,1}{:,nFeatures},'Vartype','equal');
        else 
        [tTestHypothesis,pValueTtest] = ttest2(globalMeanFeatures{1,1}{:,nFeatures},globalMeanFeatures{2,1}{:,nFeatures}, 'Vartype','unequal');
        end
        
    else 
        tTestHypothesis = NaN;
        pValueTtest = NaN;
    end
    
        allTtestHypothesis = [allTtestHypothesis, mergevars(table(tTestHypothesis, pValueTtest), [1 2], 'NewVariableName',featuresNames{nFeatures,1})];
end    

%% Save variables and export to excel
allShapiroWilkHypothesis = [table(allFilesName, 'VariableName', {'Condition'}), allShapiroWilkHypothesis];
save(fullfile(files(1).folder, 'statisticTests.mat'), 'allTtestHypothesis', 'allShapiroWilkHypothesis','featuresVarianceHypothesis')
writetable(allTtestHypothesis, fullfile(files(1).folder,'statisticalTtest.xls'),'Range','B2');
system('taskkill /F /IM EXCEL.EXE');
writetable(featuresVarianceHypothesis, fullfile(files(1).folder,'statisticalTtest.xls'),'Range','B5');
writetable(allShapiroWilkHypothesis, fullfile(files(1).folder,'statisticalShapiroWilkTest.xls'),'Range','B5');
