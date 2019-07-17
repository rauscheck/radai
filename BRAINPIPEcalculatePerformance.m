function [topPctCorrect,top3PctCorrect,outtable] = WMcalculatePerformance(pathToInputFiles, varargin)
% Read in the output of the Bayesian network and calculate performance
%
%  [topPctCorrect,top3PctCorrect,outtable] = WMcalculatePerformance(pathToInputFiles,varargin)
%
%         case 'exclude': [ acc numbers ]; case numbers to exclude
%             
%         case 'BayesOutputFN': filename, in CSV format, and it is the output of the Bayesian Network listing probabilities for all diseases for each case
%             
%         case 'CorrectDiagnosesFN': filename, in XLSX format, and has a column labeled ACC and a column labeled DIAGNOSIS
% 
%         case 'FeatureFN': full path to input file of Bayesian network.
%         If provided, it will ouput the features used for each case into
%         the output table
%
% pathToInputFiles is the root directory  
%

while ~isempty(varargin)
    switch lower(varargin{1})
        case 'exclude'
            excludeCases = varargin{2};
            
        case 'bayesoutputfn'
            BayesOutputFN = varargin{2};
            
        case 'correctdiagnosesfn'
            CorrectDiagnosesFN = varargin{2};
            
        case 'featurefn'
            FeatureFN = varargin{2};
            featuresUsed = readtable(FeatureFN);
            featuresUsed.Properties.VariableNames{1}='Acc';
            outputFeatures = 1;
            
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2)=[];
end

if ~exist('excludeCases','var')
    excludeCases = [];
end

if ~exist('pathToInputFiles','var')
    pathToInputFiles = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/'
end
cd(pathToInputFiles)

if ~exist('BayesOutputFN','var')
    BayesOutputFN = 'result_WM.csv';
end

if ~exist('FeatureFN','var')
    FeatureFN = '/Users/Andreas/Dropbox/Research/WhiteMatter-Dropbox/Thresholding/BayesInput.csv';
    outputFeatures = 0;
end

BayesOutput = readtable(BayesOutputFN);

if ~exist('CorrectDiagnosesFN','var')
    CorrectDiagnosesFN = 'correctDiagnoses-training.xlsx';
end

cd(pathToInputFiles)
correctDiagnoses = readtable(CorrectDiagnosesFN);

% initialize
totalNumCases = size(BayesOutput,1);
numTopDiagnosis = 0;
numTop3 = 0;
numExcludedCases = 0;
nonNormalCount = 0;

TopCorrect = [];
Top3Correct = [];

for subjnum=1:size(BayesOutput,1)
    curAcc = BayesOutput.ACC(subjnum);
    curDiagnosis = correctDiagnoses.DIAGNOSIS{correctDiagnoses.ACC==curAcc};
    curDiagnosis = lower(curDiagnosis);
    
    if strcmp(curDiagnosis,'normal') || sum(curAcc == excludeCases)>0 % exclude normal cases and excludeCases
        totalNumCases = totalNumCases-1;
        numExcludedCases = numExcludedCases+1;
        CorrectDiagnosis{subjnum}=upper(curDiagnosis);
        Diagnosis1{subjnum}='NA';
        Diagnosis2{subjnum}='NA';
        Diagnosis3{subjnum}='NA';
        Percent1(subjnum)=0;
        Percent2(subjnum)=0;
        Percent3(subjnum)=0;
        Acc(subjnum)=curAcc;
        TopCorrect = [TopCorrect NaN];
        Top3Correct = [Top3Correct NaN];
    else  % for all non-normal cases, calculate percent correct  
        nonNormalCount = nonNormalCount+1; % for making table at end
        
        % Get machine diagnosis
        curPercents = table2array(BayesOutput(subjnum,2:end));
        [curSortedPercents, curDiseaseInds] = sort(curPercents,'descend');
        curDiseaseInds = curDiseaseInds+1; % to index into the variable names, need to account for ACC column
        machineTop3 = {BayesOutput.Properties.VariableNames{curDiseaseInds(1:3)}};
        machineTop3 = lower(machineTop3);
        machineTop3pct = curSortedPercents(1:3);
        
        % Compare to correct diagnosis
        if strcmp(machineTop3{1},curDiagnosis)  % top diagnosis
            numTopDiagnosis = numTopDiagnosis+1;
            numTop3 = numTop3+1;
            TopCorrect = [TopCorrect 1];
            Top3Correct = [Top3Correct 1];
        elseif (strcmp(machineTop3{2},curDiagnosis) || strcmp(machineTop3{3},curDiagnosis))
            numTop3 = numTop3+1;
            TopCorrect = [TopCorrect 0];
            Top3Correct = [Top3Correct 1];
        else  % didn't get it in top 3
            TopCorrect = [TopCorrect 0];
            Top3Correct = [Top3Correct 0];
        end
        
        % Also create a table of top 3 diagnoses and percents, and if correct
        CorrectDiagnosis{subjnum}=upper(curDiagnosis);
        Diagnosis1{subjnum}=upper(machineTop3{1});
        Diagnosis2{subjnum}=upper(machineTop3{2});
        Diagnosis3{subjnum}=upper(machineTop3{3});
        Percent1(subjnum)=machineTop3pct(1);
        Percent2(subjnum)=machineTop3pct(2);
        Percent3(subjnum)=machineTop3pct(3);
        Acc(subjnum)=curAcc;
        
    end
    
end

% create the final table (excluding normals and excludeCases)
Acc=Acc';
CorrectDiagnosis=CorrectDiagnosis';
Diagnosis1=Diagnosis1';
Diagnosis2=Diagnosis2';
Diagnosis3=Diagnosis3';
Percent1=Percent1';
Percent2=Percent2';
Percent3=Percent3';
TopCorrect=TopCorrect';
Top3Correct=Top3Correct';

outtable = table(Acc,CorrectDiagnosis,Diagnosis1,Diagnosis2,Diagnosis3,Percent1,Percent2,Percent3,TopCorrect,Top3Correct);

if outputFeatures
    outtable = join(outtable,featuresUsed);  % will join based on common "Acc" column
end

outFN = 'WhAM-Diagnoses.xlsx';
if exist(outFN,'file')   % otherwise will append rather than replace, so delete the file first
    cmd=['rm ',outFN];
    system(cmd)
end
writetable(outtable,outFN,'WriteVariableNames',1)
fprintf('\n\nSaved output table to: ''WhAM-Diagnoses.xlsx'' (in current directory)\n')

% output percent correct
fprintf('Number of excluded cases, other than normals: %d\n',length(excludeCases));
fprintf('Number of excluded cases, including normals: %d\n',numExcludedCases);
fprintf('Total number of cases used: %d\n\n',totalNumCases);
topPctCorrect = numTopDiagnosis/totalNumCases
top3PctCorrect = numTop3/totalNumCases