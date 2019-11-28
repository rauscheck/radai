% Get top 1, 2, and 3 performance for AI and for humans
answersFN = 'pathToHumanAnswers.xlsx';
AIanswers = 'pathToAIAnswers.xlsx';

[AIacc1,AIacc2,AIacc3,TPR,diseaseListAI,accuracyListAI,accOrder]=WMaccuracyAI(AIanswers);  % TPR is by disease
[acc1r(1),acc2r(1),acc3r(1),TPRr1,diseaseListH]=WMaccuracyHuman(answersFN,'ResidentsCombined','WM','dd');
[acc1r(2),acc2r(2),acc3r(2),TPRr2,diseaseListH]=WMaccuracyHuman(answersFN,'ResidentsCombined','WM','ap');
[acc1r(3),acc2r(3),acc3r(3),TPRr3,diseaseListH]=WMaccuracyHuman(answersFN,'ResidentsCombined','WM','tm');
[acc1r(4),acc2r(4),acc3r(4),TPRr4,diseaseListH]=WMaccuracyHuman(answersFN,'ResidentsCombined','WM','ab');

[acc1f(1),acc2f(1),acc3f(1),TPRf1,diseaseListH]=WMaccuracyHuman(answersFN,'Fellows','WM','mc');
[acc1f(2),acc2f(2),acc3f(2),TPRf2,diseaseListH]=WMaccuracyHuman(answersFN,'Fellows','WM','ac');

[acc1c(1),acc2c(1),acc3c(1),TPRc1,diseaseListH]=WMaccuracyHuman(answersFN,'CommunityAttendings','WM','ak');
[acc1c(2),acc2c(2),acc3c(2),TPRc2,diseaseListH]=WMaccuracyHuman(answersFN,'CommunityAttendings','WM','je');

[acc1a(1),acc2a(1),acc3a(1),TPRa1,diseaseListH]=WMaccuracyHuman(answersFN,'AcademicAttendings','WM','sm');
[acc1a(2),acc2a(2),acc3a(2),TPRa2,diseaseListH]=WMaccuracyHuman(answersFN,'AcademicAttendings','WM','in');

fprintf('\n\nTOP 1')
fprintf('\nMean resident performance: %0.1f\n', mean(acc1r)*100)
fprintf('SEM resident performance: %0.1f\n', (std(acc1r)/sqrt(4))*100)
fprintf('R1: %0.1f\n',acc1r(1)*100)
fprintf('R2: %0.1f\n',acc1r(2)*100)
fprintf('R3: %0.1f\n',acc1r(3)*100)
fprintf('R4: %0.1f\n',acc1r(4)*100)
fprintf('\nMean general radiologist performance: %0.1f\n', mean(acc1c)*100)
fprintf('SEM general radiologist performance: %0.1f\n', (std(acc1c)/sqrt(2))*100)
fprintf('Gen1: %0.1f\n',acc1c(1)*100)
fprintf('Gen2: %0.1f\n',acc1c(2)*100)
fprintf('\nMean neuro fellow performance: %0.1f\n', mean(acc1f)*100)
fprintf('SEM neuro fellow performance: %0.1f\n', (std(acc1f)/sqrt(2))*100)
fprintf('F1: %0.1f\n',acc1f(1)*100)
fprintf('F2: %0.1f\n',acc1f(2)*100)
fprintf('\nMean academic performance: %0.1f\n', mean(acc1a)*100)
fprintf('SEM academic performance: %0.1f\n', (std(acc1a)/sqrt(2))*100)
fprintf('Acad1: %0.1f\n',acc1a(1)*100)
fprintf('Acad2: %0.1f\n',acc1a(2)*100)
fprintf('\nAI performance: %0.1f\n', AIacc1*100)


fprintf('\n\nTOP 2')
fprintf('\nMean resident performance: %0.1f\n', mean(acc2r)*100)
fprintf('SEM resident performance: %0.1f\n', (std(acc2r)/sqrt(4))*100)
fprintf('R1: %0.1f\n',acc2r(1)*100)
fprintf('R2: %0.1f\n',acc2r(2)*100)
fprintf('R3: %0.1f\n',acc2r(3)*100)
fprintf('R4: %0.1f\n',acc2r(4)*100)
fprintf('\nMean general radiologist performance: %0.1f\n', mean(acc2c)*100)
fprintf('SEM general radiologist performance: %0.1f\n', (std(acc2c)/sqrt(2))*100)
fprintf('Gen1: %0.1f\n',acc2c(1)*100)
fprintf('Gen2: %0.1f\n',acc2c(2)*100)
fprintf('\nMean neuro fellow performance: %0.1f\n', mean(acc2f)*100)
fprintf('SEM neuro fellow performance: %0.1f\n', (std(acc2f)/sqrt(2))*100)
fprintf('F1: %0.1f\n',acc2f(1)*100)
fprintf('F2: %0.1f\n',acc2f(2)*100)
fprintf('\nMean academic performance: %0.1f\n', mean(acc2a)*100)
fprintf('SEM academic performance: %0.1f\n', (std(acc2a)/sqrt(2))*100)
fprintf('Acad1: %0.1f\n',acc2a(1)*100)
fprintf('Acad2: %0.1f\n',acc2a(2)*100)
fprintf('\nAI performance: %0.1f\n', AIacc2*100)

fprintf('\n\nTOP 3')
fprintf('\nMean resident performance: %0.1f\n', mean(acc3r)*100)
fprintf('SEM resident performance: %0.1f\n', (std(acc3r)/sqrt(4))*100)
fprintf('R1: %0.1f\n',acc3r(1)*100)
fprintf('R2: %0.1f\n',acc3r(2)*100)
fprintf('R3: %0.1f\n',acc3r(3)*100)
fprintf('R4: %0.1f\n',acc3r(4)*100)
fprintf('\nMean general radiologist performance: %0.1f\n', mean(acc3c)*100)
fprintf('SEM general radiologist performance: %0.1f\n', (std(acc3c)/sqrt(2))*100)
fprintf('Gen1: %0.1f\n',acc3c(1)*100)
fprintf('Gen2: %0.1f\n',acc3c(2)*100)
fprintf('\nMean neuro fellow performance: %0.1f\n', mean(acc3f)*100)
fprintf('SEM neuro fellow performance: %0.1f\n', (std(acc3f)/sqrt(2))*100)
fprintf('F1: %0.1f\n',acc3f(1)*100)
fprintf('F2: %0.1f\n',acc3f(2)*100)
fprintf('\nMean academic performance: %0.1f\n', mean(acc3a)*100)
fprintf('SEM academic performance: %0.1f\n', (std(acc3a)/sqrt(2))*100)
fprintf('Acad1: %0.1f\n',acc3a(1)*100)
fprintf('Acad2: %0.1f\n',acc3a(2)*100)
fprintf('\nAI performance: %0.1f\n', AIacc3*100)

%% STATS (specifically, McNemar's test to compare each radiologist to AI system)
fprintf('\n\n----STATS USING MCNEMAR''S TEST (P VALUES COMPARING AI TO EACH RAD)----- \n\n')
allDDxPath='/Users/andreas/Documents/Work/PennRadiology/Research/Gee/WhiteMatterNetwork/HumanValidation/ValidationHumanData/AllDDx.xlsx';
rads={'resAB','resTM','resAP','resDD','genAK','genJE','felMC','felAC','acSM','acIN'};

AIDDX=readtable(allDDxPath,'Sheet','AI');
correctAnswers=readtable(allDDxPath,'Sheet','CorrectAnswers');

% create accuracy matrix for AI
numcases=size(correctAnswers,1);
accuracyListAI=zeros(numcases,3);  % first column is top1, second column is top2, third column is top3

for curCase = 1:numcases % number of cases
    if strcmpi(AIDDX.DDX1{curCase},correctAnswers.CorrectDiagnosis{curCase})
        accuracyListAI(curCase,1)=1;
        accuracyListAI(curCase,2)=1;
        accuracyListAI(curCase,3)=1;
    elseif strcmpi(AIDDX.DDX2{curCase},correctAnswers.CorrectDiagnosis{curCase})
        accuracyListAI(curCase,2)=1;
        accuracyListAI(curCase,3)=1;
    elseif strcmpi(AIDDX.DDX3{curCase},correctAnswers.CorrectDiagnosis{curCase})
        accuracyListAI(curCase,3)=1;
    end
end

for currad=1:length(rads)
    radDDX=readtable(allDDxPath,'Sheet',rads{currad});
    if sum(radDDX.ACC==AIDDX.Acc)~=length(correctAnswers.Acc) % check to make sure the accession numbers align
        error('Different order of accession numbers in radDDx vs. AI DDx suspected. Please check.')
    end
    for curCase = 1:numcases % number of cases
        if strcmpi(radDDX.DDX1{curCase},correctAnswers.CorrectDiagnosis{curCase})
            accuracyList(currad,curCase,1)=1;
            accuracyList(currad,curCase,2)=1;
            accuracyList(currad,curCase,3)=1;
        elseif strcmpi(radDDX.DDX2{curCase},correctAnswers.CorrectDiagnosis{curCase})
            accuracyList(currad,curCase,2)=1;
            accuracyList(currad,curCase,3)=1;
        elseif strcmpi(radDDX.DDX3{curCase},correctAnswers.CorrectDiagnosis{curCase})
            accuracyList(currad,curCase,3)=1;
        end
    end
end


fprintf('TOP 1 DIAGNOSIS\n')
for currad=1:length(rads)
    %currad=1;

    %[h,p,e1,e2] = testcholdout(upper(AIDDX.DDX1),upper(radDDX.DDX1),upper(correctAnswers.CorrectDiagnosis));
    [h,p,e1,e2] = testcholdout(accuracyList(currad,:,1),accuracyListAI(:,1),ones(numcases,1));
    fprintf('%s: %0.4f\n',rads{currad},p);
end

fprintf('\nTOP 2 DIAGNOSIS\n')

for currad=1:length(rads)
    %currad=1;
    %radDDX=readtable(allDDxPath,'Sheet',rads{currad});

    %numcases=length(radDDX.ACC);
    %accuracyList=zeros(length(rads),numcases,3);  % first column is top1, second column is top2, third column is top3


    [h,p,e1,e2] = testcholdout(accuracyList(currad,:,2),accuracyListAI(:,2),ones(numcases,1));
    fprintf('%s: %0.4f\n',rads{currad},p);
end


fprintf('\nTOP 3 DIAGNOSIS\n')
for currad=1:length(rads)
    %currad=1;
    %radDDX=readtable(allDDxPath,'Sheet',rads{currad});
    %numcases=length(radDDX.ACC);
    %accuracyList=zeros(numcases,3);  % first column is top1, second column is top2, third column is top3

%     for curCase = 1:numcases % number of cases
%         if strcmpi(radDDX.DDX1{curCase},correctAnswers.CorrectDiagnosis{curCase})
%             accuracyList(curCase,1)=1;
%             accuracyList(curCase,2)=1;
%             accuracyList(curCase,3)=1;
%         elseif strcmpi(radDDX.DDX2{curCase},correctAnswers.CorrectDiagnosis{curCase})
%             accuracyList(curCase,2)=1;
%             accuracyList(curCase,3)=1;
%         elseif strcmpi(radDDX.DDX3{curCase},correctAnswers.CorrectDiagnosis{curCase})
%             accuracyList(curCase,3)=1;
%         end
%     end

    [h,p,e1,e2] = testcholdout(accuracyList(currad,:,3),accuracyListAI(:,3),ones(numcases,1));
    fprintf('%s: %0.4f\n',rads{currad},p);
end



%% OUTPUT ACCURACYLIST AS CSV FILE FOR GEE MODEL
% create properly formatted table/CSV from accuracyList and accuracyListAI
% there will be 92*11 rows (11 for each radiologist/AI)
% columns will be:
%  individual: 1-11, with 11 being AI
%  experience: 1(res), 2(gen), 3(fel), 4(att), 5(AI)
%  diseasePrev: 1(rare),2(mod),3(common)
%  caseNum: 1-92 for each individual
%  disease: name of disease from ground truth (in case we want it sometime)
%accLong=reshape(accuracyList,numcases*length(rads),3);
% can do below with reshape, but this will guarantee correct dimensions
numrows=numcases*(length(rads)+1);
correctIncorrect=zeros(numrows,3);
individual=zeros(numrows,1);
experience=zeros(numrows,1);
for currad=1:(length(rads))
    curmat=accuracyList(currad,:,:);
    curmax=(currad)*numcases;
    curmin=curmax-numcases+1;
    correctIncorrect(curmin:curmax,:)=squeeze(curmat);
    individual(curmin:curmax)=currad;
end

% experience:
% individuals 1,2,3,4 are residents
% individuals 5,6 are general
% individuals 7,8 are fellows
% individuals 9,10 are attendings
% individual 11 is AI
experience(individual==1)=1;
experience(individual==2)=1;
experience(individual==3)=1;
experience(individual==4)=1;
experience(individual==5)=2;
experience(individual==6)=2;
experience(individual==7)=3;
experience(individual==8)=3;
experience(individual==9)=4;
experience(individual==10)=4;

% diseasePrevalence
diseasePrevalence = zeros(numrows,1);
diseaseText = repmat(correctAnswers.CorrectDiagnosis,length(rads)+1,1);
prevKey=readtable(allDDxPath,'Sheet','Prevalence');
for xx=1:size(prevKey,1)  % for each disease
    curDiag=prevKey.Disease{xx}; % get current diagnosis
    keyind=strcmpi(prevKey.Disease,curDiag);
    curPrev=prevKey.Prevalence{keyind};  % get prevalence of this diagnosis
    if strcmpi(curPrev,'Common')
        diseasePrevalence(strcmpi(diseaseText,curDiag))=3;
    elseif strcmpi(curPrev,'Moderate')
        diseasePrevalence(strcmpi(diseaseText,curDiag))=2;
    elseif strcmpi(curPrev,'Rare')
        diseasePrevalence(strcmpi(diseaseText,curDiag))=1;
    else
        error('Unrecognized prevalence category')
    end
end

% caseNum
caseNum=repmat([1:numcases]',length(rads)+1,1);

% add the AI results at the end
correctIncorrect(curmax+1:curmax+numcases,:)=accuracyListAI;
individual(curmax+1:curmax+numcases,:)=ones(numcases,1)*11;
experience(curmax+1:curmax+numcases,:)=ones(numcases,1)*5;

% get separate columns of 0s and 1s for resident, fellow, etc
numrows=numcases*(length(rads)+1); % +1 for AI
General=zeros(numrows,1);
Resident=zeros(numrows,1);
Fellow=zeros(numrows,1);
Attending=zeros(numrows,1);
AI=zeros(numrows,1);
for xcase=1:(numrows)
    if experience(xcase)==1
        Resident(xcase)=1;
    elseif experience(xcase)==2
        General(xcase)=1;
    elseif experience(xcase)==3
        Fellow(xcase)=1;
    elseif experience(xcase)==4
        Attending(xcase)=1;
    elseif experience(xcase)==5
        AI(xcase)=1;
    end
end

% create entire table and save out
outtable=table(caseNum,individual,experience,Resident,General,Fellow,Attending,AI,diseaseText,diseasePrevalence,correctIncorrect);
outtablePath='outpath.xlsx'; % format required for GEE
writetable(outtable,outtablePath)
