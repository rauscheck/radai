function [acc1,acc2,acc3,TPR,diseaseList,accuracyList,casesAcc] = WMaccuracyAI(answersFN,ploton)
% Calculates accuracy for top1, top2, and top3 performance for AI
%
% amr 5/13/19
%

if ~exist('ploton')
    ploton=0;
end

answerData = readtable(answersFN);

fprintf('\nCalculating accuracies for AI pipeline\n')

% make it specific for WM or BG, if asked to do so
WMBG='WMBG';  % just take all the answers

caseIdx = ~strcmpi(answerData.CorrectDiagnosis,'N/A');  % just get all the cases, but could theoretically exclude some here
cases = answerData(caseIdx,:);
% get a list of all the diseases
diseaseList = unique(upper(answerData.CorrectDiagnosis));

% exclude some random stuff that could sneak in, as well as normals (for WM)
nmind = strcmpi(diseaseList,'NORMAL') | strcmpi(diseaseList,'WM AND BG SPECIFIC PERFORMANCE') | ...
    strcmpi(diseaseList,'OVERALL PERFORMANCE') | strcmpi(diseaseList,'');
diseaseList = diseaseList(~nmind);

%cm = zeros(length(diseaseList),length(diseaseList),3);  % initialize
sens = zeros(length(diseaseList),3);
spec = zeros(length(diseaseList),3);
AUC = zeros(length(diseaseList),1);
numCases = zeros(length(diseaseList),1);

% For each disease, calculate sens, spec, AUC
for curDiseaseNum = 1:size(diseaseList,1) % for each real disease
    curDisease = diseaseList{curDiseaseNum};
    % make a new table of all times where the real diagnosis was curDisease
    foo = strcmpi(cases.CorrectDiagnosis,curDisease);  % case insensitive, index instances where correct diagnosis was current disease
    curDiseaseTable = cases(foo,:);
    curOtherDiseasesTable = cases(~foo,:);  % table of all the cases that weren't that disease (for FP and TN)
    numCases(curDiseaseNum)=size(curDiseaseTable,1);

    % for each disease, get sensitivity and specificity for that subject,
    % then average across diseases
    curmatches1 = strcmpi(curDiseaseTable.Diagnosis1,curDisease);
    curmatches2 = strcmpi(curDiseaseTable.Diagnosis2,curDisease);
    curmatches3 = strcmpi(curDiseaseTable.Diagnosis3,curDisease);

    curTP1=sum(curmatches1);
    curFN1=sum(~curmatches1);
    curTP2=sum(curmatches2 | curmatches1);
    curFN2=sum(~curmatches2 & ~curmatches1);
    curTP3=sum(curmatches3 | curmatches2 | curmatches1);
    curFN3=sum(~curmatches3 & ~curmatches2 & ~curmatches1);

    curwrongmatches1 = strcmpi(curOtherDiseasesTable.Diagnosis1,curDisease);
    curwrongmatches2 = strcmpi(curOtherDiseasesTable.Diagnosis2,curDisease);
    curwrongmatches3 = strcmpi(curOtherDiseasesTable.Diagnosis3,curDisease);

    curFP1=sum(curwrongmatches1);
    curFP2=sum(curwrongmatches1 | curwrongmatches2);
    curFP3=sum(curwrongmatches1 | curwrongmatches2 | curwrongmatches3);

    curTN1=sum(~curwrongmatches1);
    curTN2=sum(~curwrongmatches1 & ~curwrongmatches2);
    curTN3=sum(~curwrongmatches1 & ~curwrongmatches2 & ~curwrongmatches3);
    sens(curDiseaseNum,1)=curTP1/(curTP1+curFN1);  % sens = TP/(TP+FN)
    spec(curDiseaseNum,1)=curTN1/(curTN1+curFP1);  % spec = TN/(TN+FP)
    sens(curDiseaseNum,2)=curTP2/(curTP2+curFN2);  % sens = TP/(TP+FN)
    spec(curDiseaseNum,2)=curTN2/(curTN2+curFP2);  % spec = TN/(TN+FP)
    sens(curDiseaseNum,3)=curTP3/(curTP3+curFN3);  % sens = TP/(TP+FN)
    spec(curDiseaseNum,3)=curTN3/(curTN3+curFP3);  % spec = TN/(TN+FP)
    TPR(curDiseaseNum,1:5) = [0 sens(curDiseaseNum,1) sens(curDiseaseNum,2) sens(curDiseaseNum,3) 1];
    FPR(curDiseaseNum,1:5) = [0 1-spec(curDiseaseNum,1) 1-spec(curDiseaseNum,2) 1-spec(curDiseaseNum,3) 1];
    %figure; plot(FPR(curDiseaseNum,1:5),TPR(curDiseaseNum,1:5)); title(curDisease)
    %axis([0 1 0 1])
    AUC(curDiseaseNum)=trapz(FPR(curDiseaseNum,:),TPR(curDiseaseNum,:));
end

% Average TPRs and FPRs across diseases
meanTPR = mean(TPR);
meanFPR = mean(FPR);
if ploton
    figure; plot(meanFPR,meanTPR); title(sprintf('Mean ROC Across Diseases for AI pipeline'))
end

%% All above code is taken from WMAIroc.  Now to just calculate accuracies directly... note this is same as TPR except that we won't average across ...
%% diseases for this purpose.  We will just take accuracies across all cases.
numcases = size(answerData,1);
acc1=0;
acc2=0;
acc3=0;

% topcorrectInds=strcmpi(answerData.Diagnosis1,answerData.CorrectDiagnosis);
% top2Inds=strcmpi(answerData.Diagnosis2,answerData.CorrectDiagnosis);
% top3Inds=strcmpi(answerData.Diagnosis3,answerData.CorrectDiagnosis);

% going to do this via for loop in case the same diagnosis comes up twice
% (only applicable to humans, but still)
% also allows making a list (length numcases) of 1s and 0s if correct or
% incorrect-- used later for McNemar's test.  First column of matrix is for
% top1, second column is top2, third column is top3

accuracyList=zeros(numcases,3); %initialization

for curCase = 1:numcases % number of cases
    if strcmpi(answerData.Diagnosis1{curCase},answerData.CorrectDiagnosis{curCase})
        acc1=acc1+1;
        acc2=acc2+1;
        acc3=acc3+1;
        accuracyList(curCase,1)=1;
        accuracyList(curCase,2)=1;
        accuracyList(curCase,3)=1;
    elseif strcmpi(answerData.Diagnosis2{curCase},answerData.CorrectDiagnosis{curCase})
        acc2=acc2+1;
        acc3=acc3+1;
        accuracyList(curCase,2)=1;
        accuracyList(curCase,3)=1;
    elseif strcmpi(answerData.Diagnosis3{curCase},answerData.CorrectDiagnosis{curCase})
        acc3=acc3+1;
        accuracyList(curCase,3)=1;
    end
end

acc1=acc1/numcases;
acc2=acc2/numcases;
acc3=acc3/numcases;

overallAUC = trapz(meanFPR,meanTPR);

casesAcc=cases.Acc;
