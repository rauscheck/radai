% Script to complete analyze the output of the WM pipeline (WhAMpipe)
%
% 1) Thresholding of WhAMpipe output into format for Bayesian Network
% 2) Run the Bayesian Network (using python)
% 3) calculate performance comparing to correct diagnosis
% 4) Plot confusion matrix
%

desktop = 0; % requires different python command due to how it was installed; % 0 is laptop
training = 0;  % 1 for training, 0 for testing/validation, 2 for ideal cases
attendingFeatures = 0;  % 0 or 1 --> 0 means WhAMpipe features

% generally the flags below should be set to 0; only set to 1 for special analyses
CNN = 1; % if using CNN segmentation, set this flag to 1 (only applies to testing cases)
manseg = 0; % if using manual lesion segmentations, set this flag to 1; also need to choose training or validation under training flag above

%% Exclude the following cases, for technical purposes (registration problems, etc)
if training==1 && ~attendingFeatures  %training
    excludeCases = []; %[26639772; 5569269; 26744538; 27028898; 28923987; 15654929; 8881731; 15035975; 2687494];
elseif training==0 && ~attendingFeatures % testing
    excludeCases = []; %[14284233; 14380031; 14789179; 15040132; 15372988; 15537513; 23993000; 26954347; 27223618; 27239765; 27647635; 6498571; 6835454]; ...
    %14096571; 14745970; 15035975; 15484399; 15654929; 15879532; 26753078; 26846939; 2687494; 26906282; 27951991; 28227356; 28923987; 3043357; 6204101; 7725020; 8743450; 8765883];  % the second line is all ones where T1 to FLAIR registration was fixed manually
else
    excludeCases = [];
end



%% Thresholding of WhAMpipe output into format for Bayesian Network
% To do any adjustment of thresholds, open up WMthresholding and edit at the top
if CNN==1
    pathToInputFiles = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/Thresholding/CNN';
    FeatureFile = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/WMfeatures-CNN.csv';
    BayesOutFile = 'result_WM_CNN.csv';
    CorrectDiagnosesFN = 'correctDiagnoses-allCases.xlsx';
    
elseif manseg==1
    if training==1  % training
        pathToInputFiles = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/Thresholding/ManSeg/Training';
        FeatureFile = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/WMfeatures-manseg.csv';
        BayesOutFile = 'result_WM_manseg-training.csv';
        CorrectDiagnosesFN = 'correctDiagnoses-allCases.xlsx';
        
    elseif training==0  % validation cases
        pathToInputFiles = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/Thresholding/ManSeg/Validation';
        FeatureFile = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/WMfeatures-manseg.csv';
        BayesOutFile = 'result_WM_manseg-validation.csv';
        CorrectDiagnosesFN = 'correctDiagnoses-allCases.xlsx';
        
    else
        error('Invalid parameter for "training" variable');
    end
    
else
    
    if training==1  % training
        if attendingFeatures
            FeatureFile = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/WMfeatures-attending_training.csv';
            CorrectDiagnosesFN = 'correctDiagnoses-attending_training.xlsx';
            BayesOutFile = 'result_WM_attending_training.csv';
        else
            pathToInputFiles = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/Thresholding/Training';
            FeatureFile = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/WMfeatures-training.csv';
            CorrectDiagnosesFN = 'correctDiagnoses-training.xlsx';
            BayesOutFile = 'result_WM_training.csv';
        end
        
        
    elseif training==2  % ideal cases
        pathToInputFiles = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/Thresholding/Ideal';
        FeatureFile = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/WMidealCaseFeatures.csv';
        CorrectDiagnosesFN = 'correctDiagnoses-idealWMcases.csv';
        
        
    elseif training==0  % validation
        if ~attendingFeatures
            pathToInputFiles = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/Thresholding/Validation';
            FeatureFile = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/WMfeatures-validation.csv';
            BayesOutFile = 'result_WM_validation.csv';
            CorrectDiagnosesFN = 'correctDiagnoses-validation.xlsx';
        elseif attendingFeatures
            FeatureFile = '/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/WMfeatures-attending_validation.csv';
            CorrectDiagnosesFN = 'correctDiagnoses-attending_validation.xlsx';
            BayesOutFile = 'result_WM_attending_validation.csv';
        end
    else
        error('Invalid parameter for "training" variable');
    end
    
    
end

%% Do the Thresholding
if ~attendingFeatures % this part is to create thresholded features, so don't run for attending features
    if training==1 || training==0  % don't run for ideal cases (which already have "thresholded features")
        addpath('/Users/Andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/Thresholding');
        threshStats = WMthresholding(pathToInputFiles,FeatureFile);
    end
end

%% Run the Bayesian Network (using python)
cd('/Users/andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/NaiveBayes/');

% For laptop
if ~desktop
    system(sprintf(['/anaconda/envs/python2/bin/python2.7 SimpleNaiveBayes.py WMBayesProbs_simple.xlsx "v6" %s %s'],FeatureFile,BayesOutFile));
    %system(sprintf(['/anaconda/envs/python2/bin/python2.7 SimpleNaiveBayes.py WMBayesProbs_simple.xlsx "v5simple" %s %s'],FeatureFile,BayesOutFile));
end
% For desktop
if desktop
    %system(sprintf(['/Users/andreas/anaconda2/bin/python SimpleNaiveBayes.py WMBayesProbs_simple.xlsx "v5simple" %s %s'],FeatureFile,BayesOutFile));
    system(sprintf(['/Users/andreas/anaconda2/bin/python SimpleNaiveBayes.py WMBayesProbs_simple.xlsx "v6" %s %s'],FeatureFile,BayesOutFile));
end

%!/anaconda/envs/python2/bin/python2.7 SimpleNaiveBayes.py WMBayesProbs_simple.xlsx "v5simple" BayesInput.csv result_WM_training.csv

% Run for ideal cases scenario:
%!/anaconda/envs/python2/bin/python2.7 SimpleNaiveBayes.py WMBayesProbs_simple.xlsx "v5simple" WMidealCaseFeatures.csv result_WM_training.csv


%% Calculate performance comparing to correct diagnosis
[topPctCorrect,top3PctCorrect,outtable] = WMcalculatePerformance(pwd,'exclude',excludeCases,'FeatureFN',FeatureFile,'BayesOutputFN',BayesOutFile,'CorrectDiagnosesFN',CorrectDiagnosesFN); % for normal processing

% Run for ideal cases scenario:
%[topPctCorrect,top3PctCorrect,outtable] = WMcalculatePerformance(pwd,'BayesOutputFN','result_WM_training.csv','CorrectDiagnosesFN','correctDiagnoses-idealWMcases.csv');
%% Plot a confusion matrix for all diseases-- think about taking out "normals"

% get a list of all the diseases
diseaseList = unique(outtable.CorrectDiagnosis);

% exclude normals (for white matter only)
nmind = strcmpi(diseaseList,'NORMAL');
diseaseList = diseaseList(~nmind);

cm = zeros(length(diseaseList),length(diseaseList));  % initialize
% confusionTable = outtable;
% confusionTable{strcmp(confusionTable.CorrectDiagnosis(:)=='Normal'}='';

% For each disease, fill in the confusion matrix
for curDiseaseNum = 1:size(diseaseList,1) % for each real disease
    curDisease = diseaseList{curDiseaseNum};
    % make a new table of all times where the real diagnosis was curDisease
    foo = strcmpi(outtable.CorrectDiagnosis,curDisease);  % case insensitive, index instances where correct diagnosis was current disease
    curDiseaseTable = outtable(foo,:);
    
    % for each disease, get percent of time the correct disease was called the predicted disease
    for curPredNum = 1:size(diseaseList,1)
        curPred = diseaseList{curPredNum};
        curmatches = strcmpi(curDiseaseTable.Diagnosis1,curPred);
        curPct = sum(curmatches)/size(curDiseaseTable,1);  % percent of times predicted disease matched
        cm(curPredNum,curDiseaseNum)=curPct;  % rows are predicted disease, columns are real disease (columns add up to 1)
    end
end

% Plot the confusion matrix, and output the list of diseases in order
plotConfusionMatrix(cm,diseaseList);  % optional third argument to leave axis labels off
