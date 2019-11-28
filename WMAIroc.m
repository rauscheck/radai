function [sens,spec,TPR,FPR,overallAUC,partAUC,targetsLong,confidencesLong] = WMAIroc(answersFN,ploton)
% Plots overall AUC curve for AI algorithm from a particular excel sheet.
% Similar to WMroc, but with input excel sheet format adapted.
%
% ploton is a flag (0 or 1) to determine whether to plot a figure of ROC
%
% excludes normal cases-- not sure what happens if you try to include them...
%
%
% amr 12/26/18
%

if ~exist('ploton')
    ploton=0;
end

answerData = readtable(answersFN);

fprintf('\nPrinting ROC for AI pipeline\n')

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

        % try an alternative way of calculating ROCs by getting confidence scores and "targets" as in matlab function roc
    numTargetCases=length(curmatches1);
    numNonTargets=length(curwrongmatches1);
    numTotalCases=numTargetCases+numNonTargets;
    targets(curDiseaseNum,1:numTargetCases)=1;  % vector of numbers with targets being 1, others as 0s
    targets(curDiseaseNum,numTargetCases+1:numTargetCases+numNonTargets)=0;
    % confidences (0=target not present, 1=ddx3, 2=ddx2, 3=ddx1, so that highest confidence is for DDX1)
    %confidences = zeros(numTotalCases,1);
    confidences(curDiseaseNum,curmatches1)=3;
    confidences(curDiseaseNum,curmatches2)=2;
    confidences(curDiseaseNum,curmatches3)=1;
    wrong1=logical([zeros(numTargetCases,1);curwrongmatches1]);  % make it length of total number cases; vector of where the disease was called the top 1 diagnosis but wrongly
    wrong2=logical([zeros(numTargetCases,1);curwrongmatches2]);  % make it length of total number cases
    wrong3=logical([zeros(numTargetCases,1);curwrongmatches3]);  % make it length of total number cases
    if sum(wrong1)>0, confidences(curDiseaseNum,wrong1)=3; end   % only do if there are actual indices
    if sum(wrong2)>0, confidences(curDiseaseNum,wrong2)=2; end
    if sum(wrong3)>0, confidences(curDiseaseNum,wrong3)=1; end
end

% Average TPRs and FPRs across diseases
meanTPR = mean(TPR);
meanFPR = mean(FPR);
if ploton
    figure; plot(meanFPR,meanTPR); title(sprintf('Mean ROC Across Diseases for AI pipeline'))
end
overallAUC = trapz(meanFPR,meanTPR);
partAUC = trapz(meanFPR(1:4),meanTPR(1:4))/meanFPR(4); % calculate AUC up to the maximal specifici

% Targets and confidence scores -- this produces almost identical results to averaging across diseases, as I have done above.
numrows=size(targets,1)*size(targets,2);
targetsLong=reshape(targets,numrows,1);
confidencesLong=reshape(confidences,numrows,1);
