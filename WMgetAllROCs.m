% Compare different confusion matrices

averages=0; % if 1, then just plot average non-parametric ROC curve. if 0, then will do parametric ROC curves

answersFN = 'pathToHumanAnswers.xlsx';

close all
figure; hold on
linw =3;
marksz=1;

rescolor =[0.2 0.2 0.2]; %[0.1 0.3 0.3]; %[0.5 0.2 0.2];
felcolor =[0.5 0.5 0.5];
comcolor = [0.3 0.3 0.3];
accolor =[0.6 0.6 0.6]; %[0.5 0.3 0.3]; %[0.6 0.6 0.6];
AIcolor= [0.1 0.6 0.1];


AIanswers = 'pathToAIAnswers.xlsx';
[sens,spec,TPR,FPR,overallAUC_AI,partAUC_AI,targetsAI,confidencesAI] = WMAIroc(AIanswers);
if averages==1
    ai=plot(mean(FPR),mean(TPR),'-ko','LineWidth',linw+4,...
        'MarkerSize',marksz,'Color',AIcolor);
end



if averages==0

    % parametric ROC for AI
    x=confidencesAI(logical(targetsAI))+0.1;  % add 1 because they need to be postivie scores
    y=confidencesAI(logical(~targetsAI))+0.1;
    models='normal-normal'; %'weibull-weibull'; %'normal-normal'; %'gamma-gamma'; %
    paucrange=[0 0.3];  % for partial AUC, set to [] if you do not want partial AUC
    %[rc aucparam aucpartial]=paramroc(y,x,models,paucrange,1);
    t=[0:0.01:1];
    showNonParametric=0;  % this needs to be 1 if you want empirical AUC values (aucemp)
    nbootstrap=100;
    confInterval=[]; %t; %[]; %t
    provideConfidenceIntervals=1;
    lnwdth=3;
    %[hand rc aucparam aucpartial aucemp CIareaparam CIareapartial CIareaemp CIparamroc]=paramroc(x,y,models,fprrange,np,bootssams,atwhich)
    %[hAI rc aucparam aucpartial aucemp CIareaparam CIareapartial CIareaemp CIparamroc]=paramroc(y,x,models,[],0,100,t);
    %[hAI rc aucparam aucpartial aucemp CIareaparam CIareapartial CIareaemp CIparamroc]=paramroc(y,x,models,[],1,100,t);

    % show confidence intervals
    if provideConfidenceIntervals
        [hAI rc aucparam aucpart_AI aucemp_AI CIareaparam_AI CIareapartial_AI CIareaemp_AI CIparamroc_AI]=paramroc(y,x,models,paucrange,showNonParametric,nbootstrap,confInterval);
    else
        [hAI rc aucparam aucpart_AI aucemp_AI]=paramroc(y,x,models,[],showNonParametric,nbootstrap,confInterval);
    end

    hAI.LineWidth=lnwdth;

    hold on;
    [sens,spec,TPR,FPR,overallAUC_r1,partAUC_r1,targets_r1,confidences_r1] = WMroc(answersFN,'ResidentsCombined','WM','dd');
    %r1=plot(mean(FPR),mean(TPR),'-ro','LineWidth',linw,...
    %    'MarkerSize',marksz);
    [sens,spec,TPR,FPR,overallAUC_r2,partAUC_r2,targets_r2,confidences_r2] = WMroc(answersFN,'ResidentsCombined','WM','ap');
    %r2=plot(mean(FPR),mean(TPR),'-ro','LineWidth',linw,...
    %    'MarkerSize',marksz);
    [sens,spec,TPR,FPR,overallAUC_r3,partAUC_r3,targets_r3,confidences_r3] = WMroc(answersFN,'ResidentsCombined','WM','tm');
    %r3=plot(mean(FPR),mean(TPR),'-ro','LineWidth',linw,...
    %    'MarkerSize',marksz);
    [sens,spec,TPR,FPR,overallAUC_r4,partAUC_r4,targets_r4,confidences_r4] = WMroc(answersFN,'ResidentsCombined','WM','ab');
    %r4=plot(mean(FPR),mean(TPR),'-ro','LineWidth',linw,...
    %    'MarkerSize',marksz);

    % resident ROC (parametric), with option to do partial auc and bootstraps
    resTargets = [targets_r1;targets_r2;targets_r3;targets_r4];
    resConfidences = [confidences_r1;confidences_r2;confidences_r3;confidences_r4];
    x=resConfidences(logical(resTargets))+0.1;  % add 1 because they need to be postivie scores
    y=resConfidences(logical(~resTargets))+0.1;
    if provideConfidenceIntervals
        [hres rc aucparam aucpart_res aucemp_res CIareaparam_res CIareapartial_res CIareaemp_res CIparamroc_res]=paramroc(y,x,models,paucrange,showNonParametric,nbootstrap,confInterval);
    else
        [hres rc aucparam aucpart_res aucemp_res]=paramroc(y,x,models,[],showNonParametric,nbootstrap,confInterval);
    end
    hres.LineWidth=lnwdth;


    [sens,spec,TPR,FPR,overallAUC_f1,partAUC_f1,targets_f1,confidences_f1] = WMroc(answersFN,'Fellows','WM','mc');
    %f1=plot(mean(FPR),mean(TPR),'--go','LineWidth',linw,...
    %    'MarkerSize',marksz);
    [sens,spec,TPR,FPR,overallAUC_f2,partAUC_f2,targets_f2,confidences_f2] = WMroc(answersFN,'Fellows','WM','ac');
    %f2=plot(mean(FPR),mean(TPR),'--go','LineWidth',linw,...
    %    'MarkerSize',marksz);

    % fellow ROC (parametric), with option to do partial auc and bootstraps
    hold on
    felTargets = [targets_f1;targets_f2];
    felConfidences = [confidences_f1;confidences_f2];
    x=felConfidences(logical(felTargets))+0.1;  % add 1 because they need to be postivie scores
    y=felConfidences(logical(~felTargets))+0.1;
    if provideConfidenceIntervals
        [hfel rc aucparam aucpart_fel aucemp_fel CIareaparam_fel CIareapartial_fel CIareaemp_fel CIparamroc_fel]=paramroc(y,x,models,paucrange,showNonParametric,nbootstrap,confInterval);
    else
        [hfel rc aucparam aucpart_fel aucemp_fel]=paramroc(y,x,models,[],showNonParametric,nbootstrap,confInterval);
    end
    hfel.LineWidth=lnwdth;

    [sens,spec,TPR,FPR,overallAUC_c1,partAUC_c1,targets_c1,confidences_c1] = WMroc(answersFN,'CommunityAttendings','WM','ak');
    %c1=plot(mean(FPR),mean(TPR),':bo','LineWidth',linw,...
    %    'MarkerSize',marksz);
    [sens,spec,TPR,FPR,overallAUC_c2,partAUC_c2,targets_c2,confidences_c2] = WMroc(answersFN,'CommunityAttendings','WM','je');
    %c2=plot(mean(FPR),mean(TPR),':bo','LineWidth',linw,...
    %    'MarkerSize',marksz);

    % community parametric ROC
    hold on
    comTargets = [targets_c1;targets_c2];
    comConfidences = [confidences_c1;confidences_c2];
    x=comConfidences(logical(comTargets))+0.1;  % add 1 because they need to be postivie scores
    y=comConfidences(logical(~comTargets))+0.1;
    if provideConfidenceIntervals
        [hc rc aucparam aucpart_gen aucemp_gen CIareaparam_gen CIareapartial_gen CIareaemp_gen CIparamroc_gen]=paramroc(y,x,models,paucrange,showNonParametric,nbootstrap,confInterval);
    else
        [hc rc aucparam aucpart_gen aucemp_gen]=paramroc(y,x,models,[],showNonParametric,nbootstrap,confInterval);
    end
    hc.LineWidth=lnwdth;

    % attendings
    [sens,spec,TPR,FPR,overallAUC_a1,partAUC_a1,targets_a1,confidences_a1] = WMroc(answersFN,'AcademicAttendings','WM','sm');
    %a1=plot(mean(FPR),mean(TPR),'-.mo','LineWidth',linw,...
    %    'MarkerSize',marksz);
    [sens,spec,TPR,FPR,overallAUC_a2,partAUC_a2,targets_a2,confidences_a2] = WMroc(answersFN,'AcademicAttendings','WM','in');
    %a2=plot(mean(FPR),mean(TPR),'-.mo','LineWidth',linw,...
    %    'MarkerSize',marksz);

    % attending parametric ROC
    hold on
    attTargets = [targets_a1;targets_a2];
    attConfidences = [confidences_a1;confidences_a2];
    x=attConfidences(logical(attTargets))+0.1;  % add 1 because they need to be postivie scores
    y=attConfidences(logical(~attTargets))+0.1;
    if provideConfidenceIntervals
        [hatt rc aucparam aucpart_att aucemp_att CIareaparam_att CIareapartial_att CIareaemp_att CIparamroc_att]=paramroc(y,x,models,paucrange,showNonParametric,nbootstrap,confInterval);
    else
        [hatt rc aucparam aucpart_att aucemp_att]=paramroc(y,x,models,[],showNonParametric,nbootstrap,confInterval);
    end
    hatt.LineWidth=lnwdth;

    fprintf('Resident 1 AUC: %0.3f\n',overallAUC_r1)
    fprintf('Resident 2 AUC: %0.3f\n',overallAUC_r2)
    fprintf('Resident 3 AUC: %0.3f\n',overallAUC_r3)
    fprintf('Resident 4 AUC: %0.3f\n',overallAUC_r4)
    fprintf('Mean: %0.3f\n\n',mean([overallAUC_r1 overallAUC_r2 overallAUC_r3 overallAUC_r4]))

    fprintf('Community Attending 1 AUC: %0.3f\n',overallAUC_c1)
    fprintf('Community Attending 2 AUC: %0.3f\n',overallAUC_c2)
    fprintf('Mean: %0.3f\n\n',mean([overallAUC_c1 overallAUC_c2]))

    fprintf('Fellow 1 AUC: %0.3f\n',overallAUC_f1)
    fprintf('Fellow 2 AUC: %0.3f\n',overallAUC_f2)
    fprintf('Mean: %0.3f\n\n',mean([overallAUC_f1 overallAUC_f2]))

    fprintf('Academic Attending 1 AUC: %0.3f\n',overallAUC_a1)
    fprintf('Academic Attending 2 AUC: %0.3f\n',overallAUC_a2)
    fprintf('Mean: %0.3f\n\n',mean([overallAUC_a1 overallAUC_a2]))

    fprintf('AI pipeline AUC: %0.3f\n\n',overallAUC_AI)

    % parametric partial AUCs based on paucrange
    fprintf('PARTIAL AUCs\n')
    fprintf('Resident partial AUC: %0.3f [%0.3f %0.3f]\n',aucpart_res/paucrange(2),CIareapartial_res./paucrange(2))
    fprintf('General rad partial AUC: %0.3f [%0.3f %0.3f]\n',aucpart_gen/paucrange(2),CIareapartial_gen./paucrange(2))
    fprintf('Fellow partial AUC: %0.3f [%0.3f %0.3f]\n',aucpart_fel/paucrange(2),CIareapartial_fel./paucrange(2))
    fprintf('Attending partial AUC: %0.3f [%0.3f %0.3f]\n',aucpart_att/paucrange(2),CIareapartial_att./paucrange(2))
    fprintf('AI pipeline partial AUC: %0.3f [%0.3f %0.3f]\n',aucpart_AI/paucrange(2),CIareapartial_AI./paucrange(2))

elseif averages==1
    [sens,spec,TPRr1,FPRr1,overallAUC_r1,partAUC_r1] = WMroc(answersFN,'ResidentsCombined','WM','dd');
    [sens,spec,TPRr2,FPRr2,overallAUC_r2,partAUC_r2] = WMroc(answersFN,'ResidentsCombined','WM','ap');
    [sens,spec,TPRr3,FPRr3,overallAUC_r3,partAUC_r3] = WMroc(answersFN,'ResidentsCombined','WM','tm');
    [sens,spec,TPRr4,FPRr4,overallAUC_r4,partAUC_r4] = WMroc(answersFN,'ResidentsCombined','WM','ab');
    meanTPRr=mean([mean(TPRr1);mean(TPRr2);mean(TPRr3);mean(TPRr4)]);
    meanFPRr=mean([mean(FPRr1);mean(FPRr2);mean(FPRr3);mean(FPRr4)]);
    r1=plot(meanFPRr,meanTPRr,'-ro','LineWidth',linw,'MarkerSize',marksz);
    resSTD=std([mean(TPRr1);mean(TPRr2);mean(TPRr3);mean(TPRr4)]);
    %r1=shadedErrorBar(meanFPRr,meanTPRr,resSTD,'lineProps',{'-ro','Color',rescolor,'LineWidth',linw,'MarkerSize',marksz});

    [sens,spec,TPRf1,FPRf1,overallAUC_f1,partAUC_f1] = WMroc(answersFN,'Fellows','WM','mc');
    [sens,spec,TPRf2,FPRf2,overallAUC_f2,partAUC_f2] = WMroc(answersFN,'Fellows','WM','ac');
    meanTPRf=mean([mean(TPRf1);mean(TPRf2)]);
    meanFPRf=mean([mean(FPRf1);mean(FPRf2)]);
    f1=plot(meanFPRf,meanTPRf,'--go','LineWidth',linw,'MarkerSize',marksz);
    felSTD=std([mean(TPRf1);mean(TPRf2)]);
    %f1=shadedErrorBar(meanFPRf,meanTPRf,felSTD,'lineProps',{'--go','Color',felcolor,'LineWidth',linw,'MarkerSize',marksz});

    [sens,spec,TPRc1,FPRc1,overallAUC_c1,partAUC_c1] = WMroc(answersFN,'CommunityAttendings','WM','ak');
    [sens,spec,TPRc2,FPRc2,overallAUC_c2,partAUC_c2] = WMroc(answersFN,'CommunityAttendings','WM','je');
    meanTPRc=mean([mean(TPRc1);mean(TPRc2)]);
    meanFPRc=mean([mean(FPRc1);mean(FPRc2)]);
    c1=plot(meanFPRc,meanTPRc,':bo','LineWidth',linw,'MarkerSize',marksz);
    cSTD=std([mean(TPRc1);mean(TPRc2)]);
    %c1=shadedErrorBar(meanFPRc,meanTPRc,cSTD,'lineProps',{':bo','Color',comcolor,'LineWidth',linw,'MarkerSize',marksz});

    [sens,spec,TPRa1,FPRa1,overallAUC_a1,partAUC_a1] = WMroc(answersFN,'AcademicAttendings','WM','sm');
    [sens,spec,TPRa2,FPRa2,overallAUC_a2,partAUC_a2] = WMroc(answersFN,'AcademicAttendings','WM','in');
    meanTPRa=mean([mean(TPRa1);mean(TPRa2)]);
    meanFPRa=mean([mean(FPRa1);mean(FPRa2)]);
    a1=plot(meanFPRa,meanTPRa,'-.mo','LineWidth',linw,'MarkerSize',marksz);
    aSTD=std([mean(TPRa1);mean(TPRa2)]);
    %a1=shadedErrorBar(meanFPRa,meanTPRa,aSTD,'lineProps',{'-.mo','Color',accolor,'LineWidth',linw,'MarkerSize',marksz});
end


% Plot formatting
xlabel('1-Specificity','FontSize',16)
ylabel('Sensitivity','FontSize',16)
%xlabel('False Positive Rate','FontSize',16)
%ylabel('True Positive Rate','FontSize',16)
ax = gca;
ax.FontSize=30;
ax.FontName='Arial';
ax.FontWeight='bold';
f = gcf;

% Set better color style of lines

if ~averages
        %hres.Color = rescolor;
        %hfel.Color = [0 0 0.6]; %felcolor;
        %hc.Color = [0.2 0.6 0.6]; %comcolor;
        hatt.Color = [0.2 0.2 0.2];
        hAI.Color = [0 0.6 0];
        hAI.LineWidth=lnwdth+4;

        % style
        %hres.LineStyle='--'
    %     r1.Color = rescolor; r2.Color = rescolor; r3.Color = rescolor; r4.Color = rescolor;
    %     f1.Color = felcolor; f2.Color = felcolor;
    %     c1.Color = comcolor; c2.Color = comcolor;
    %     a1.Color = accolor; a2.Color = accolor;
else
    r1.Color = rescolor;
    f1.Color = felcolor;
    c1.Color = comcolor;
    a1.Color = accolor;
    ai.Color = AIcolor;
end

% 0 to 1 line
hold on
l=plot([0 1],[0 1],':ko','LineWidth',2);
l.Color=[0.7 0.7 0.7];

if averages

    % turn markers off
    if ~averages
        r1.Marker='none';r2.Marker='none';r3.Marker='none';r4.Marker='none';
        f1.Marker='none';f2.Marker='none';
        c1.Marker='none';c2.Marker='none';
        a1.Marker='none';a2.Marker='none';
    else
        r1.Marker='none';
        f1.Marker='none';
        c1.Marker='none';
        a1.Marker='none';
    end
    ai.Marker='none';

    [lgd,icons]=legend([r1 c1 f1 a1 ai],...
        '   Radiology Residents (AUC 0.728)','   General Radiologists (AUC 0.709)',...
        '   Neuroradiology Fellows (AUC 0.845)','   Academic Neuroradiologists (AUC 0.902)',...
        '   AI System (AUC 0.912)');
    lgd.Location='southeast';
    lgd.FontSize=16;
    lgd.FontName='Arial';
    lgd.Color=[0.93 0.93 0.93];
    lgd.Position=[0.2798 0.1739 0.6 0.2];
    %lgd.Box='off';
else
    %allPlots=get(gca,'Children');
    % new empirical AUC measuresement (slightly different)
    [lgd,icons]=legend([hres hc hfel hatt hAI],...
                '   Radiology Residents (AUC 0.731)','   General Radiologists (AUC 0.718)',...
        '   Neuroradiology Fellows (AUC 0.849)','   Academic Neuroradiologists (AUC 0.904)',...
        '   AI algorithm (AUC 0.920)');
    %legend(allPlots([5, 6, 3, 2, 4, 1]), {'AI', 'Academic', 'Fellows', 'General', 'Residents', 'Chance'});

end

% legend appearance
lgd.FontSize=14;
lgd.FontName='Arial';
lgd.Color=[0.93 0.93 0.93];
lgd.Position=[0.4 0.1739 0.52 0.2];

% resize figure
pos = get(f,'position');
set(f,'position',[pos(1:2)*2 pos(3:4)*2])

% other formatting
if ~averages
    title('')
end
