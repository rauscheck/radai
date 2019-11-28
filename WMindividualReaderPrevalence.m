alldata='pathToGEEformat.xlsx';
d=readtable(alldata);

% number of top3 correct for each category

% initialization for total numbers
nCor_allradCommon=0;
nCor_allradMod=0;
nCor_allradRare=0;
nTot_allradCommon=0;
nTot_allradMod=0;
nTot_allradRare=0;


%% loop over individuals
% loop over prevalence
for prev=1:3  % common =3, rare =1

    for radNumber=1:10 %11 total readers, with 11 being AI
        idx1=d.individual==radNumber;
        idx2=d.diseasePrevalence==prev; % common =3, rare =1

        curIntersect=(idx1&idx2);  % intersection of idx1 (individual) with idx2 (prevalence)
        nCor=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
        nTot=size(d.correctIncorrect_3(curIntersect),1);  % total number

        %fprintf('Rad number %0.0f prevalence %d: Correct %d Total %d Perc %0.1f \n',radNumber, prev,nCor,nTot,nCor/nTot*100)
        fprintf('%0.1f ',nCor/nTot*100)
    end
end
