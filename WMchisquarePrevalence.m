alldata='pathToData.xlsx';
d=readtable(alldata);

% number of top3 correct for each category

% initialization for total numbers
nCor_allradCommon=0;
nCor_allradMod=0;
nCor_allradRare=0;
nTot_allradCommon=0;
nTot_allradMod=0;
nTot_allradRare=0;


%% residents
idx1=d.Resident==1;

% common
idx2=d.diseasePrevalence==3; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_resCommon=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_resCommon=size(d.correctIncorrect_3(curIntersect),1);  % total number

nCor_allradCommon=nCor_allradCommon+nCor_resCommon;
nTot_allradCommon=nTot_allradCommon+nTot_resCommon;

% moderate
idx2=d.diseasePrevalence==2; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_resMod=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_resMod=size(d.correctIncorrect_3(curIntersect),1);  % total number

nCor_allradMod=nCor_allradMod+nCor_resMod;
nTot_allradMod=nTot_allradMod+nTot_resMod;

% rare
idx2=d.diseasePrevalence==1; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_resRare=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_resRare=size(d.correctIncorrect_3(curIntersect),1);  % total number

nCor_allradRare=nCor_allradRare+nCor_resRare;
nTot_allradRare=nTot_allradRare+nTot_resRare;

% Chi Square Comparing common to rare
n1 = nCor_resCommon; N1 = nTot_resCommon;
n2 = nCor_resRare; N2 = nTot_resRare;
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2)
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected)
p_residents = 1 - chi2cdf(chi2stat,1)


%% General Radiologists
idx1=d.General==1;

% common
idx2=d.diseasePrevalence==3; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Common=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Common=size(d.correctIncorrect_3(curIntersect),1);  % total number

% moderate
idx2=d.diseasePrevalence==2; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Mod=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Mod=size(d.correctIncorrect_3(curIntersect),1);  % total number

% rare
idx2=d.diseasePrevalence==1; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Rare=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Rare=size(d.correctIncorrect_3(curIntersect),1);  % total number

% Chi Square Comparing common to rare
n1 = nCor_Common; N1 = nTot_Common;
n2 = nCor_Rare; N2 = nTot_Rare;
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2)
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected)
p_general = 1 - chi2cdf(chi2stat,1)

% add for total tally
nCor_allradCommon=nCor_allradCommon+nCor_Common;
nTot_allradCommon=nTot_allradCommon+nTot_Common;

nCor_allradMod=nCor_allradMod+nCor_Mod;
nTot_allradMod=nTot_allradMod+nTot_Mod;

nCor_allradRare=nCor_allradRare+nCor_Rare;
nTot_allradRare=nTot_allradRare+nTot_Rare;


%% Fellows
idx1=d.Fellow==1;

% common
idx2=d.diseasePrevalence==3; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Common=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Common=size(d.correctIncorrect_3(curIntersect),1);  % total number

% moderate
idx2=d.diseasePrevalence==2; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Mod=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Mod=size(d.correctIncorrect_3(curIntersect),1);  % total number

% rare
idx2=d.diseasePrevalence==1; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Rare=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Rare=size(d.correctIncorrect_3(curIntersect),1);  % total number

% Chi Square Comparing common to rare
n1 = nCor_Common; N1 = nTot_Common;
n2 = nCor_Rare; N2 = nTot_Rare;
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2)
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected)
p_fellow = 1 - chi2cdf(chi2stat,1)

% add to totals for final tally
nCor_allradCommon=nCor_allradCommon+nCor_Common;
nTot_allradCommon=nTot_allradCommon+nTot_Common;

nCor_allradMod=nCor_allradMod+nCor_Mod;
nTot_allradMod=nTot_allradMod+nTot_Mod;

nCor_allradRare=nCor_allradRare+nCor_Rare;
nTot_allradRare=nTot_allradRare+nTot_Rare;


%% Academic Attending
idx1=d.Attending==1;

% common
idx2=d.diseasePrevalence==3; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Common=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Common=size(d.correctIncorrect_3(curIntersect),1);  % total number

% moderate
idx2=d.diseasePrevalence==2; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Mod=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Mod=size(d.correctIncorrect_3(curIntersect),1);  % total number

% rare
idx2=d.diseasePrevalence==1; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Rare=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Rare=size(d.correctIncorrect_3(curIntersect),1);  % total number

% Chi Square Comparing common to rare
n1 = nCor_Common; N1 = nTot_Common;
n2 = nCor_Rare; N2 = nTot_Rare;
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2)
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected)
p_academic = 1 - chi2cdf(chi2stat,1)

% add to totals for final tally
nCor_allradCommon=nCor_allradCommon+nCor_Common;
nTot_allradCommon=nTot_allradCommon+nTot_Common;

nCor_allradMod=nCor_allradMod+nCor_Mod;
nTot_allradMod=nTot_allradMod+nTot_Mod;

nCor_allradRare=nCor_allradRare+nCor_Rare;
nTot_allradRare=nTot_allradRare+nTot_Rare;


%% AI
idx1=d.AI==1;

% common
idx2=d.diseasePrevalence==3; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Common=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Common=size(d.correctIncorrect_3(curIntersect),1);  % total number

% moderate
idx2=d.diseasePrevalence==2; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Mod=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Mod=size(d.correctIncorrect_3(curIntersect),1);  % total number

% rare
idx2=d.diseasePrevalence==1; % common =3, rare =1
curIntersect=(idx1&idx2);  % intersection of idx1 (specialty) with idx2 (prevalence)
nCor_Rare=sum(d.correctIncorrect_3(curIntersect));  % get out only the correct ones
nTot_Rare=size(d.correctIncorrect_3(curIntersect),1);  % total number

% Chi Square Comparing common to rare
n1 = nCor_Common; N1 = nTot_Common;
n2 = nCor_Rare; N2 = nTot_Rare;
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2)
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
chi2stat = sum((observed-expected).^2 ./ expected)
p_AI = 1 - chi2cdf(chi2stat,1)
