function [sens,spec,TPR,FPR,overallAUC,partAUC,targetsLong,confidencesLong] = WMroc(answersFN,sheetname,WMBG,subjInitials,ploton)
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
% excludes normal cases-- not sure what happens if you try to include them...
%
% 6/10/19:  added outputs targetsLong and confidencesLong
%  targetsLong is array of 1s and 0s where 1s are TPs or FPs, and 0s are TNs and FNs
%  confidencesLong is array of same length with "scores" provided by radiologist, with 3 being DDx1, 1 being DDx3, and 0 being not in DDx (4 confidence levels)
%
% 11/21/19:  added output partAUC, which is a partial AUC calculated based
% on the maximal specificity (x-value) in the dataset, and normalized by
% this, so that AUC value is not extrapolated from that point to [1,1]
%
% amr 12/26/18; modified 6/10/19 and 11/21/19
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

fprintf('\nPrinting ROC for %s with diseases %s\n\n',sheetname,WMBG)

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

confidences=zeros(length(diseaseList),size(cases,1));
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
    wrong1=logical([zeros(numTargetCases,1);curwrongmatches1]);  % make it length of total number cases
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
figure; plot(meanFPR,meanTPR); title(sprintf('Mean Across Diseases %s %s',sheetname,upper(subjInitials)))
end
overallAUC = trapz(meanFPR,meanTPR);
partAUC = trapz(meanFPR(1:4),meanTPR(1:4))/meanFPR(4); % calculate AUC up to the maximal specificity and normalize by that

% Targets and confidence scores -- this produces almost identical results to averaging across diseases, as I have done above.
numrows=size(targets,1)*size(targets,2);
targetsLong=reshape(targets,numrows,1);
confidencesLong=reshape(confidences,numrows,1);

% this allows you to use plotroc function
[tpr,fpr,thresholds] = roc(targetsLong',confidencesLong'./3);
if ploton
    figure; plotroc(targetsLong',confidencesLong')
    x=confidencesLong(logical(targetsLong));
    y=confidencesLong(logical(~targetsLong));
    models='normal-normal'; %'gamma-gamma';
    paucrange=[0 0.3];  % for partial AUC
    [rc aucparam aucpartial]=paramroc(y,x,models,paucrange,1);
    t=[0:0.01:1];
    [rc aucparam aucpartial aucemp CIareaparam CIareapartial CIareaemp CIparamroc]=paramroc(y,x,models,paucrange,1,100,t)
end
% % Also try plotting using plotroc function
% numDiseases = length(diseaseList);
%
% numberTotalCases = length(
% targets=zeros(length(diseaseList),
% for foo=1:length(diseaseList)
%     curDz=diseaseList{foo};
%
% end

%plotroc(targets, outputs)
