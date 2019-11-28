function [acc1,acc2,acc3,TPR,diseaseList,accuracyList]=WMaccuracyHuman(answersFN,sheetname,WMBG,subjInitials,ploton)
% Plots overall AUC curves for subject participants from a particular excel sheet
%
% must specify one of following sheet names within excel file:
%       sheetname = 'AcademicAttendings';
%       sheetname = 'ResidentsCombined';
%       sheetname = 'Fellows';
%       sheetname = 'CommunityAttendings';
%
% WMBG must be 'WM' or 'BG' (or you can try 'WMBG' and hope for the best)
% subjInitials are the 2 initials of the subject in question
%
% OUTPUTS:
% acc1,2,3 are just percentages calculated across all cases
% TPR is the true positive rate for each disease separately, specified in
% diseaseList
% accuracyList is a matrix size numCases x 3, with 1s and 0s for correct/incorrect of
% cases in order, for top1 (first column), top2 (2nd column), and top3 (3rd column)
%
% excludes normal cases-- not sure what happens if you try to include them...
%
% amr 12/26/18, modified 6/5/2019
%


if ~exist('sheetname','var')
    sheetname = 'AcademicAttendings';  % other choices ResidentsCombined,Fellows,CommunityAttendings,AcademicAttendings
    %sheetname = 'ResidentsCombined';
    %sheetname = 'Fellows';
    %sheetname = 'CommunityAttendings';
end

if ~exist('WMBG','var')
    WMBG = 'WM';  % can also be 'BG' or 'WMBG'/'BGWM' (for both)
end

if ~exist('ploton','var')
    ploton=0;
end

if ~exist('subjInitials','var')
    subjInitials = 'in';
end

answerData = readtable(answersFN,'sheet',sheetname);


% make it specific for WM or BG, if asked to do so
WMcaseIdx = strcmpi(answerData.BGTop,'WM');
BGcaseIdx = strcmpi(answerData.WMTop,'BG');
WMBGcaseIdx = WMcaseIdx | BGcaseIdx;

if strcmpi(WMBG,'WM')
    cases = answerData(WMcaseIdx,:);
elseif strcmpi(WMBG,'BG')
    cases = answerData(BGcaseIdx,:);
elseif strcmpi(WMBG,'WMBG') || strcmpi(WMBG,'BGWM')
    cases = answerData(WMBGcaseIdx,:);
else
    error('Inappropriate choice of WMBG variable')
end

% make it specific to the subject initials, if asked to do so
if ~strcmpi(subjInitials,'none')
    subjIdx = strcmpi(cases.Subj,subjInitials);
    cases = cases(subjIdx,:);
end

% get a list of all the diseases
diseaseList = unique(upper(cases.CorrectAnswer));

% exclude some random stuff that could sneak in, as well as normals (for WM)
nmind = strcmpi(diseaseList,'NORMAL') | strcmpi(diseaseList,'WM AND BG SPECIFIC PERFORMANCE') | ...
    strcmpi(diseaseList,'OVERALL PERFORMANCE') | strcmpi(diseaseList,'');
diseaseList = diseaseList(~nmind);

cm = zeros(length(diseaseList),length(diseaseList),3);  % initialize
sens = zeros(length(diseaseList),3);
spec = zeros(length(diseaseList),3);
AUC = zeros(length(diseaseList),1);
numCases = zeros(length(diseaseList),1);

% For each disease, calculate sens, spec, AUC
for curDiseaseNum = 1:size(diseaseList,1) % for each real disease
    curDisease = diseaseList{curDiseaseNum};
    % make a new table of all times where the real diagnosis was curDisease
    foo = strcmpi(cases.CorrectAnswer,curDisease);  % case insensitive, index instances where correct diagnosis was current disease
    curDiseaseTable = cases(foo,:);
    curOtherDiseasesTable = cases(~foo,:);  % table of all the cases that weren't that disease (for FP and TN)
    numCases(curDiseaseNum)=size(curDiseaseTable,1);

    % for each disease, get sensitivity and specificity for that subject,
    % then average across diseases
    curmatches1 = strcmpi(curDiseaseTable.WM1,curDisease) | strcmpi(curDiseaseTable.BG1,curDisease);
    curmatches2 = strcmpi(curDiseaseTable.WM2,curDisease) | strcmpi(curDiseaseTable.BG2,curDisease);
    curmatches3 = strcmpi(curDiseaseTable.WM3,curDisease) | strcmpi(curDiseaseTable.BG3,curDisease);

    curTP1=sum(curmatches1);
    curFN1=sum(~curmatches1);
    curTP2=sum(curmatches2 | curmatches1);
    curFN2=sum(~curmatches2 & ~curmatches1);
    curTP3=sum(curmatches3 | curmatches2 | curmatches1);
    curFN3=sum(~curmatches3 & ~curmatches2 & ~curmatches1);

    curwrongmatches1 = strcmpi(curOtherDiseasesTable.WM1,curDisease) | strcmpi(curOtherDiseasesTable.BG1,curDisease);
    curwrongmatches2 = strcmpi(curOtherDiseasesTable.WM2,curDisease) | strcmpi(curOtherDiseasesTable.BG2,curDisease);
    curwrongmatches3 = strcmpi(curOtherDiseasesTable.WM3,curDisease) | strcmpi(curOtherDiseasesTable.BG3,curDisease);

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
figure; plot(meanFPR,meanTPR); title(sprintf('Mean Across Diseases %s %s',sheetname,upper(subjInitials)))
end
overallAUC = trapz(meanFPR,meanTPR);

%% All above code is taken from WMroc.  Now to just calculate accuracies directly... note this is same as TPR except that we won't average across ...
%% diseases for this purpose.  We will just take accuracies across all cases.
numcases = size(cases,1);
acc1=0;
acc2=0;
acc3=0;

% going to do this via for loop in case the person answered the same
% diagnosis twice.  This will ensure it gets counted correctly
% also allows making a list (length numcases) of 1s and 0s if correct or
% incorrect-- used later for McNemar's test.  First column of matrix is for
% top1, second column is top2, third column is top3
accuracyList=zeros(numcases,3);

for curCase = 1:numcases % number of cases
    if strcmpi(cases.WM1{curCase},cases.CorrectAnswer{curCase})
        acc1=acc1+1;
        acc2=acc2+1;
        acc3=acc3+1;
        accuracyList(curCase,1)=1;
        accuracyList(curCase,2)=1;
        accuracyList(curCase,3)=1;
    elseif strcmpi(cases.WM2{curCase},cases.CorrectAnswer{curCase})
        acc2=acc2+1;
        acc3=acc3+1;
        accuracyList(curCase,2)=1;
        accuracyList(curCase,3)=1;
    elseif strcmpi(cases.WM3{curCase},cases.CorrectAnswer{curCase})
        acc3=acc3+1;
        accuracyList(curCase,3)=1;
    end
end

acc1=acc1/numcases;
acc2=acc2/numcases;
acc3=acc3/numcases;

% Also try plotting using plotroc function
%plotroc(targets, outputs)

%casesAcc=cases.Acc;
