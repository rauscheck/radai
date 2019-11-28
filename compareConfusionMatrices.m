% Compare different confusion matrices

answersFN = 'pathToData.xlsx';
c=[0 0.8];

cmRes = plotSubjectConfusionMatrices(answersFN,'ResidentsCombined','WM');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-ResidentsCombined.png','-dpng')

cmFel = plotSubjectConfusionMatrices(answersFN,'Fellows','WM');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-FellowsCombined.png','-dpng')

cmCommunity = plotSubjectConfusionMatrices(answersFN,'CommunityAttendings','WM');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-CommunityCombined.png','-dpng')

cmAcademic = plotSubjectConfusionMatrices(answersFN,'AcademicAttendings','WM');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-AcademicCombined.png','-dpng')

cmRes1 = plotSubjectConfusionMatrices(answersFN,'ResidentsCombined','WM','ap');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-R1.png','-dpng')

cmRes2 = plotSubjectConfusionMatrices(answersFN,'ResidentsCombined','WM','dd');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-R2.png','-dpng')

cmRes3 = plotSubjectConfusionMatrices(answersFN,'ResidentsCombined','WM','ab');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-R3.png','-dpng')

cmRes4 = plotSubjectConfusionMatrices(answersFN,'ResidentsCombined','WM','tm');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-R4.png','-dpng')

% make a perfect confusion matrix
cmPerfect = plotSubjectConfusionMatrices(answersFN,'Perfect','WM');
caxis(c)

ResToFel = corr2(cmRes,cmFel)
ResToCommunity = corr2(cmRes,cmCommunity)
ResToAcademic = corr2(cmRes,cmAcademic)
FelToAcademic = corr2(cmFel,cmAcademic)
FelToCommunity = corr2(cmFel,cmCommunity)
communityToAcademic = corr2(cmCommunity,cmAcademic)

% cm comes from running WMBayesProcessing, or can use "perfect confusion
% matrix" also
if ~exist('cm')
    cm = cmPerfect;
end
machineToRes = corr2(cm./max(max(cm)),cmRes./max(max(cmRes)))  %corr2(cm,cmRes)
machineToFel = corr2(cm./max(max(cm)),cmFel./max(max(cmFel)))  %corr2(cm,cmFel)
machineToCommunity = corr2(cm./max(max(cm)),cmCommunity./max(max(cmCommunity)))
machineToAcademic = corr2(cm./max(max(cm)),cmAcademic./max(max(cmAcademic)))

resToFel = corr2(cmRes./max(max(cmRes)),cmFel./max(max(cmFel)))
resToCommunity = corr2(cmRes./max(max(cmRes)),cmCommunity./max(max(cmCommunity)))
resToAcademic = corr2(cmRes./max(max(cmRes)),cmAcademic./max(max(cmAcademic)))


% Plot some more individual subject confusion matrices as well for intracategory
% correlations
cmCom1 = plotSubjectConfusionMatrices(answersFN,'CommunityAttendings','WM','ak');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-Community1.png','-dpng')

cmCom2 = plotSubjectConfusionMatrices(answersFN,'CommunityAttendings','WM','je');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-Community2.png','-dpng')

cmFel1 = plotSubjectConfusionMatrices(answersFN,'Fellows','WM','ac');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-Fellow1.png','-dpng')

cmFel2 = plotSubjectConfusionMatrices(answersFN,'Fellows','WM','mc');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-Fellow2.png','-dpng')

cmAcad1 = plotSubjectConfusionMatrices(answersFN,'AcademicAttendings','WM','sm');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-Academic1.png','-dpng')

cmAcad2 = plotSubjectConfusionMatrices(answersFN,'AcademicAttendings','WM','in');
caxis(c)
print('Confusion_ColorBar01/ConfusionMatrix-Academic2.png','-dpng')

FelToFel = corr2(cmFel1,cmFel2)
%ComToCom = corr2(cmCom1,cmCom2)
AcadToAcad = corr2(cmAcad1,cmAcad2)

r1tor2 = corr2(cmRes1,cmRes2);
r1tor3 = corr2(cmRes1,cmRes3);
r1tor4 = corr2(cmRes1,cmRes4);
r2tor3 = corr2(cmRes2,cmRes3);
r3tor4 = corr2(cmRes3,cmRes4);
RtoRavg = mean([r1tor2 r1tor3 r1tor4 r2tor3 r3tor4])

cmtor1 = corr2(cmRes1,cm)
cmtor2 = corr2(cmRes2,cm)
cmtor3 = corr2(cmRes3,cm)
cmtor4 = corr2(cmRes4,cm)
cmtoRavg = mean([cmtor1 cmtor2 cmtor3 cmtor4])

cmtoResCombo = corr2(cmRes,cm)

%% New statistical approach.  Turn matrix into a vector, take out all the identity values (so left only with errors)
% Compare Attendings to each other
numDiseases = size(cmAcad1,1);
cmPerfect_long=reshape(cmPerfect,[numDiseases*numDiseases 1]);

cmAcad1_long=reshape(cmAcad1,[numDiseases*numDiseases 1]);
cmAcad1_long(cmPerfect_long==1)=NaN;  % get rid of correct responses, only left with errors
cmAcad1_long=cmAcad1_long(~isnan(cmAcad1_long));

cmAcad2_long=reshape(cmAcad2,[numDiseases*numDiseases 1]);
cmAcad2_long(cmPerfect_long==1)=NaN;  % get rid of correct responses, only left with errors
cmAcad2_long=cmAcad2_long(~isnan(cmAcad2_long));

cmRes1_long=reshape(cmRes1,[numDiseases*numDiseases 1]);
cmRes1_long(cmPerfect_long==1)=NaN;
cmRes1_long=cmRes1_long(~isnan(cmRes1_long));

cmRes2_long=reshape(cmRes2,[numDiseases*numDiseases 1]);
cmRes2_long(cmPerfect_long==1)=NaN;
cmRes2_long=cmRes2_long(~isnan(cmRes2_long));

cmRes3_long=reshape(cmRes3,[numDiseases*numDiseases 1]);
cmRes3_long(cmPerfect_long==1)=NaN;
cmRes3_long=cmRes3_long(~isnan(cmRes3_long));

cmRes4_long=reshape(cmRes4,[numDiseases*numDiseases 1]);
cmRes4_long(cmPerfect_long==1)=NaN;
cmRes4_long=cmRes4_long(~isnan(cmRes4_long));

cmFel1_long=reshape(cmFel1,[numDiseases*numDiseases 1]);
cmFel1_long(cmPerfect_long==1)=NaN;
cmFel1_long=cmFel1_long(~isnan(cmFel1_long));

cmFel2_long=reshape(cmFel2,[numDiseases*numDiseases 1]);
cmFel2_long(cmPerfect_long==1)=NaN;
cmFel2_long=cmFel2_long(~isnan(cmFel2_long));

cmAI_long=reshape(cm,[numDiseases*numDiseases 1]);
cmAI_long(cmPerfect_long==1)=NaN;
cmAI_long=cmAI_long(~isnan(cmAI_long));

% get rho for various correlations
[r_acadacad,PVAL] = corr(cmAcad1_long,cmAcad2_long)
[r_felfel,PVAL] = corr(cmFel1_long,cmFel2_long)

[r_acad1fel1,PVAL] = corr(cmAcad1_long,cmFel1_long);
[r_acad1fel2,PVAL] = corr(cmAcad1_long,cmFel2_long);
[r_acad2fel1,PVAL] = corr(cmAcad2_long,cmFel1_long);
[r_acad2fel2,PVAL] = corr(cmAcad2_long,cmFel2_long);
r_avg_acadfel=mean([r_acad1fel1 r_acad1fel2 r_acad2fel1 r_acad2fel2])


[r_acad1,PVAL] = corr(cmAcad1_long,cmAI_long);
[r_acad2,PVAL] = corr(cmAcad2_long,cmAI_long);
r_avg_acadAI=mean([r_acad1 r_acad2])

[r_fel1,PVAL] = corr(cmFel1_long,cmAI_long);
[r_fel2,PVAL] = corr(cmFel2_long,cmAI_long);
r_avg_felAI=mean([r_fel1 r_fel2])

[r_res1,PVAL] = corr(cmRes1_long,cmAI_long);
[r_res2,PVAL] = corr(cmRes2_long,cmAI_long);
[r_res3,PVAL] = corr(cmRes3_long,cmAI_long);
[r_res4,PVAL] = corr(cmRes4_long,cmAI_long);
r_avg_resAI=mean([r_res1 r_res2 r_res3 r_res4])
