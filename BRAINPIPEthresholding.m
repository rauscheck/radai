function stats = BRAINPIPEthresholding(pathToInputFiles,BayesFile,saveOutputs)
% Function to convert the output of the BRain Automated INference pipeline
% (BRAINpipe) to input for the Bayesian network (BRAINnet)
%
%     stats = BRAINPIPEthresholding(pathToInputFiles,BayesFile,saveOutputs)
%           default path:  ~/default/path
%           BayesFile is the filename (full path) to the output file (input to subsequent Bayes Net)
%           saveOutputs is an optional argument (default 1) about whether to save the outputs of the function to disk
%
%  Please note that input files must have standardized names and be located within the same path:
%       Output/[acc]_IndividualLesionMeasures.csv (for each subject)-- has information on each lesion
%       Output/SummaryLesionMeasures.csv -- each row is a subject (first column specifies acc)
%       clinicalInfo.xlsx must have all the clinical features on each subject
%


%%
if ~exist('saveOutputs','var')  % optional input to choose whether to save the output of the function (or just return the stats)
    saveOutputs=1;
end

%% SET THRESHOLDS
clusterThresh = 10; % general cluster threshold for lesions to count as real lesions

T1high_volthresh = 10;  % total volume of high T1 signal within lesions (low threshold)
T1low_volPctThresh = 0.2;  % percent of lesion voxels with low T1 signal (higher threshold)
T2thresh = [-2 1];   % decreased, normal, and increased
GREvolthresh = 150;  % volume of abnormal susceptibility within lesions
enhanceThresh = 0.40; % if greater than this on any lesion, 'yes', otherwise 'none', together with volume
enhanceSizeThresh = 50; % we only take relatively larger lesions seriously for enhancement ratio
enhanceVolIndeterminate = [10 1000];  % in this range, call it indeterminate
DWIvolthresh = [100 100000];  % don't believe incredibly large values
ADCvolthresh = [100 100000];

CCthreshPctRange = [0.016 0.033];  % if more than x% of lesion volume is in CC, count it
CCthreshVolRange = [500 1000];  % if more than 1000 voxels in CC, count it
cortexThreshPctRange = [0.1 0.7]; 
periventThresh = 5; % distance to ventricles
periventVolumePctThreshRange = [0.2 0.4]; 
juxtaThresh = 2; % distance of lesions judged to be juxtacortical
juxtaPctThreshRange = [0.49 0.50]; % in between, it's indeterminate  what percentage of lesions within the distance threshold above to meet juxtacortical metric 'yes'
symmetryThresh = 0.55; % bigger than this is called asymmetric
sizeThresh = [1000 15000]; % for max size, bounds for small, medium, large
numberThreshRange = [3 6]; % less than or equal to, we'll just call it "single" lesion; greater, we call it "multiple" (could look at lesion sizes or distribution...); will only count lesions larger than 20mm^3 as lesions
antTempThresh = 10; % number of voxels found in anterior temporal lobe lesions, needs to be bilateral --> 'yes' or 'no'
antTempAsymmetryThresh = 0.3; % i.e. if 30% of voxels are only in 1 hemisphere, doesn't count for ant temp
masseffectNumberNormThresh = [0.04 0.119];  % this is mass effect ratio divided by number of lesions
frontalThresh = 60;
parietalThresh = 50;
temporalThresh = 50;
occipitalThresh = 20;

ventThreshRange = [31000 41000];   %in between is indeterminate
extentThreshRange = [10000 50000];  % direct volume of lesion burden (rather than percentage)

enhanceVolRatioThresh = [0.1 0.12]; % if less than first, then "small" ratio, if more than second, then "large" ratio


%% DO ALL PROCESSING WITHIN A SINGLE DIRECTORY, read in files
if ~exist('pathToInputFiles','var')
    pathToInputFiles = '~/default/path/'
end
cd(pathToInputFiles)

if ~exist('Output/SummaryLesionMeasures.csv','file')
    error('No file named SummaryLesionMeasures.csv found. Note you need a directory named Output that contains BRAINpipe output')
end

summaryInput = readtable('Output/SummaryLesionMeasures.csv');
clinicalInfo = readtable('clinicalInfo.xlsx');


%% INITIALIZE
% set all the output cell arrays
acc = [];  % accession numbers

% clinical features
age = [];  % non-character array
chronicity = [];
gender = [];
immunocompromised = [];
prodrome = [];

% signal features
diffusion = [];
enhancement = [];
FLAIR = [];
susceptibility = [];
T1 = [];
T2 = [];
enhancementRatio = [];

% spatial features
corpusCallosum = [];
cortex = [];
lobar_distribution = [];
mass_effect = [];
number = [];
periventricular = [];
lesionSize = [];
juxtacortical = [];
symmetry = [];
anterior_temporal = [];
ventVol = [];
lesionExtent = [];

% keep track of a few others for stats purposes
stats.CCvol = [];
stats.CCpct = [];
stats.AvgDistanceToGM = [];
stats.MaxDistanceToGM = [];
stats.MinDistanceToGM = [];
stats.NumJuxtacortLesions = [];
stats.PctJuxtacortLesions = [];
stats.AvgDistanceToVents = [];
stats.MaxDistanceToVents = [];
stats.MinDistanceToVents = [];
stats.NumPeriventLesions = [];
stats.PctPeriventLesions = [];
stats.PeriventricularVolume = [];
stats.PeriventricularVolumePct = [];
stats.MaxLesionSize = [];
stats.MinLesionSize = [];
stats.MeanLesionSize = [];
stats.MedianLesionSize = [];
stats.LesionSize85pct = [];
stats.NumLesions = [];
stats.masseffect = [];
stats.masseffectNumberNorm = [];
stats.LobarPercents = [];
stats.CortexVol = [];
stats.CortexVolPct = [];
stats.GRElowest = [];
stats.antTempTot = [];
stats.antTempR = [];
stats.antTempL = [];
stats.antTempSym = [];
stats.enhancingVol = [];
stats.VentVol = [];
stats.TotLesionVol = [];
stats.TotLesionPct = [];
stats.WMvol = [];
stats.enhanceVolRatio = [];


%% THRESHOLDING: LOOP THROUGH EACH PATIENT WITHIN SUMMARY INPUT FILE
% for each patient, calculate the new output measures with appropriate thresholding
% then start building a new table for output
for IDcount = 1:size(summaryInput,1)
    curACC=table2array(summaryInput(IDcount,'ID'));
    fprintf('Processing subject: %d \n',curACC)
    
    acc=[acc;curACC];
    
    %% CLINICAL FEATURES-- these just get read in from the clinicalFeatures table
    curage=(clinicalInfo.AGE(clinicalInfo.ACC==curACC));
    if curage<=40
        curagetext="young";
    elseif curage>=61
        curagetext="old";
    else
        curagetext="adult";
    end
    age = [age; curagetext];
    
    curchron = (clinicalInfo.CHRONICITY(clinicalInfo.ACC==curACC));
    if ~(strcmp(curchron,'acute') || strcmp(curchron,'chronic'))  %% only allow these outputs
        curchron = 'NA';
    end
    chronicity = [chronicity;curchron];
    
    curgender = (clinicalInfo.GENDER(clinicalInfo.ACC==curACC));
    if ~(strcmp(curgender,'M') || strcmp(curgender,'F'))  %% only allow these outputs
        curgender = 'NA';
    end
    gender = [gender;curgender];
    
    curimmune = (clinicalInfo.IMMUNOCOMPROMISED(clinicalInfo.ACC==curACC));
    if ~(strcmp(curimmune,'yes') || strcmp(curimmune,'no'))  %% only allow these outputs
        curimmune = 'NA';
    end
    immunocompromised = [immunocompromised;curimmune];
    
    curprodrome = (clinicalInfo.PRODROME(clinicalInfo.ACC==curACC));
    if ~(strcmp(curprodrome,'yes') || strcmp(curprodrome,'no'))  %% only allow these outputs
        curprodrome = 'NA';
    end
    prodrome = [prodrome;curprodrome];
    
    %% SIGNAL FEATURES    
    
    % overall relative T2 signal within lesions compared to background
    
    curT2 = summaryInput.Relative_T2_Signal(summaryInput.ID==curACC);
    if curT2==-1
        T2 = [T2 ; "NA"];  % may not have a T2
    else
        if curT2<=T2thresh(1)
            T2 = [T2 ; "decreased"];
        elseif curT2>=T2thresh(2)
            T2 = [T2 ; "increased"];
        else
            T2 = [T2 ; "normal"];
        end
    end
    
    FLAIR = [FLAIR ; "increased"]; % could get rid of this, but keeping just in case it's useful someday
    
    % lesion-by-lesion kinds of features
    %       load in the patient's individual lesions
    clear indivLesions
    indivFN = ['Output/' num2str(curACC) '_IndividualLesionMeasures.csv'];
    if ~exist(indivFN,'file')
        error('No IndividualLesionMeasures.csv file for %s',num2str(curACC))
    end
    indivLesions = readtable(indivFN);
    
    % T1 by number of voxels overlapping with lesions (T1 high takes precedence over T1 low)
    % We take the top 50% size lesions (smaller lesions have unreliable T1 signal), and find percentage of voxels
    % with high or low T1 within each lesion.  If lesion has any high T1 or >30% of lesion has low T1, then name it as such
    
    allLesions = indivLesions.Size_of_Lesion_mm3;
    topHalf = ceil(length(allLesions)/2);  % half of lesions
    topHalfLesionsT1high = indivLesions.T1high_vol(1:topHalf);
    topHalfLesionsT1low = indivLesions.T1low_vol(1:topHalf);
    topHalfLesionSizes = allLesions(1:topHalf);
    
    T1high = topHalfLesionsT1high >= T1high_volthresh;  % use a small volume threshold for T1high (even a small amount of high T1 signal counts)
    T1low = topHalfLesionsT1low./topHalfLesionSizes >= T1low_volPctThresh;  % lesion must have a good amount of low T1 signal to be called as such
    
    if sum(T1high)>0
        T1=[T1 ; "increased"];
    elseif sum(T1low)>0
        T1=[T1 ; "decreased"];
    else
        T1=[T1 ; "normal"];
    end
    
    
    % GRE by number of voxels in GRE U-net overlapping with lesions
    GREidx = indivLesions.GREvol >= GREvolthresh;
    if sum(GREidx)>0
        susceptibility=[susceptibility ; "yes"];
    elseif max(indivLesions.GREvol)==-1
        susceptibility=[susceptibility ; "NA"];
    else
        susceptibility=[susceptibility ; "no"];
    end
    
    
    % Enhancement by number of "enhancing" voxels overlapping with lesions
    enhanceIdx = indivLesions.Enhancement >= enhanceThresh & indivLesions.Size_of_Lesion_mm3 >= enhanceSizeThresh;
    curMaxEnhanceVol = max(indivLesions.EnhancingVolume);
    maxEnhanceRatio=max(indivLesions.Enhancement);
    if maxEnhanceRatio==-1
        enhancement=[enhancement ; "NA"];
    else
        if sum(enhanceIdx)>0 || curMaxEnhanceVol>enhanceVolIndeterminate(2)
            enhancement=[enhancement ; "yes"];
        elseif sum(enhanceIdx)>0 || curMaxEnhanceVol>enhanceVolIndeterminate(1)
            enhancement=[enhancement ; "NA"];  % indeterminate
        else
            enhancement=[enhancement ; "none"];  % definitively no enhancement
        end
    end
    stats.enhancingVol = [stats.enhancingVol curMaxEnhanceVol];
    
    % Ratio of enhancing volume to lesion size for top half size lesions
    enhLesions = find(indivLesions.EnhancingVolume>enhanceVolIndeterminate(1));
    if isempty(enhLesions)  % exclude cases with only nonenhancing lesions; include any with at least indeterminate enhancement
        curEnhanceVolRatio = -1;
        enhancementRatio=[enhancementRatio ; "NA"];
    else
        allCurRatios = indivLesions.EnhancingVolume(enhLesions)./indivLesions.Size_of_Lesion_mm3(enhLesions); % get ratio of all enhancing lesions
        meaningfulCurRatio = allCurRatios(1:ceil(end/3)); % only take top third largest lesions
        curEnhanceVolRatio = median(meaningfulCurRatio);
        if curEnhanceVolRatio>=enhanceVolRatioThresh(2)
            enhancementRatio=[enhancementRatio ; "large"];
        elseif curEnhanceVolRatio<=enhanceVolRatioThresh(1)
            enhancementRatio=[enhancementRatio ; "small"];
        else
            enhancementRatio=[enhancementRatio ; "NA"];
        end
    end
    stats.enhanceVolRatio = [stats.enhanceVolRatio curEnhanceVolRatio];
    
    
    % Diffusion
    DWIlesionIdx=indivLesions.DWIrestrictVol>DWIvolthresh(1) & indivLesions.DWIrestrictVol<DWIvolthresh(2);
    ADClesionIdx=indivLesions.ADCrestrictVol>ADCvolthresh(1) & indivLesions.ADCrestrictVol<ADCvolthresh(2);
    if min(indivLesions.Relative_DWI_Signal)==-1 && min(indivLesions.Relative_ADC_Signal)==-1  % neither exist, missing sequence
        diffusion = [diffusion ; "NA"];
    else
        if sum(DWIlesionIdx)>0  % at least 1 lesion met DWI criteria
            if sum(DWIlesionIdx & ADClesionIdx)>0  % same lesion fulfills both criteria
                diffusion = [diffusion ; "restricted"];
            elseif min(indivLesions.Relative_ADC_Signal)==-1  % if ADC doesn"t exist, but DWI fulfills criteria, still ok
                diffusion = [diffusion ; "restricted"];
            else
                diffusion = [diffusion ; "normal"];  % it was t2 shine-through
            end
        elseif sum(ADClesionIdx)>0  % at least 1 ADC lesion met criteria, but we already know DWI didn't
            if min(indivLesions.Relative_DWI_Signal)==-1  % DWI doesn't exist, so we'll just say we don't know, because ADC by itself is not sensitive/specific enough at our threshold
                diffusion = [diffusion ; "NA"];
            else  % ADC met criteria, but DWI exists and did NOT meet criteria-- probably artifact
                diffusion = [diffusion ; "normal"];
            end
        else  % if neither DWI or ADC meet criteria
            diffusion = [diffusion ; "normal"];
        end
    end
    
    
    %% SPATIAL FEATURES
    numLesions = summaryInput.Num_of_Lesions((summaryInput.ID==curACC));
    
    curCCvol = summaryInput.Total_Volume_CorpusCallosum_mm3(summaryInput.ID==curACC);
    curTotLesionVol = summaryInput.Total_Size_of_Lesions_mm3(summaryInput.ID==curACC);
    curCCpct= curCCvol/curTotLesionVol;
    stats.CCvol = [stats.CCvol curCCvol];
    stats.CCpct = [stats.CCpct curCCpct];
    
    % Corpus Callosum involvement (by percentage or volume)
    if curCCpct>=CCthreshPctRange(2) || curCCvol>=CCthreshVolRange(2)  % use either percent or volume
        corpusCallosum = [corpusCallosum ; "yes"];
    elseif curCCpct<=CCthreshPctRange(1) || curCCvol<=CCthreshVolRange(1)
        corpusCallosum = [corpusCallosum ; "no"];
    else
        corpusCallosum = [corpusCallosum ; "NA"];
    end
    
    
    curCortexVol = summaryInput.Total_Volume_GM_mm3(summaryInput.ID==curACC);
    curCortexVolPct = curCortexVol/summaryInput.Total_Size_of_Lesions_mm3(summaryInput.ID==curACC);
    stats.CortexVol = [stats.CortexVol curCortexVol];
    stats.CortexVolPct = [stats.CortexVolPct curCortexVolPct];
    
    % Cortical gray matter involvement (by percent)
    if curCortexVolPct>=cortexThreshPctRange(2)
        cortex = [cortex ; "yes"];
    elseif curCortexVolPct<=cortexThreshPctRange(1)
        cortex = [cortex ; "no"];
    else
        cortex = [cortex ; "NA"];
    end
    
    % do periventricular and juxtacortical on a per lesion basis (but also keep track of average, max, min)
    % many ways to calculate this; for stats purposes we keep track of
    % multiple measures
    clear foo
    foo = indivLesions.Dist_To_Ventricles_Voxels;
    stats.AvgDistanceToVents = [stats.AvgDistanceToVents mean(foo)];
    stats.MaxDistanceToVents = [stats.MaxDistanceToVents max(foo)];
    stats.MinDistanceToVents = [stats.MinDistanceToVents min(foo)];
    periventLesions = foo(foo<periventThresh);  % distance threshold
    numPeriventLesions = length(periventLesions);
    pctPeriventLesions = numPeriventLesions/numLesions;
    stats.NumPeriventLesions = [stats.NumPeriventLesions numPeriventLesions];
    stats.PctPeriventLesions = [stats.PctPeriventLesions pctPeriventLesions];
    
    
    % For thresholding we are doing periventricular on a volume basis (if a lesion is close to
    % ventricles, then count the volume of that lesion, sum them up, then
    % divide by total lesion volume to get percent)
    periventInds = foo<periventThresh;
    periventInds_real = indivLesions.Size_of_Lesion_mm3>clusterThresh;  % only take lesions greater than cluster threshold
    periventInds = periventInds.*periventInds_real;
    periventVolume = sum(indivLesions.Size_of_Lesion_mm3.*periventInds);
    periventVolumePct = periventVolume/summaryInput.Total_Size_of_Lesions_mm3(summaryInput.ID==curACC);
    stats.PeriventricularVolume = [stats.PeriventricularVolume periventVolume];
    stats.PeriventricularVolumePct = [stats.PeriventricularVolumePct periventVolumePct];
    
    if periventVolumePct>=periventVolumePctThreshRange(2)
        periventricular = [periventricular ; "yes"];
    elseif periventVolumePct<=periventVolumePctThreshRange(1)
        periventricular = [periventricular ; "no"];
    else
        periventricular = [periventricular ; "NA"];
    end
    
    clear foo
    foo = indivLesions.Dist_To_CorticalGM_Voxels;
    stats.AvgDistanceToGM = [stats.AvgDistanceToGM mean(foo)];
    stats.MaxDistanceToGM = [stats.MaxDistanceToGM max(foo)];
    stats.MinDistanceToGM = [stats.MinDistanceToGM min(foo)];
    juxtacortLesions = foo(foo<juxtaThresh);  % distance threshold
    numJuxtacortLesions = length(juxtacortLesions);
    pctJuxtacortLesions = numJuxtacortLesions/numLesions;
    stats.NumJuxtacortLesions = [stats.NumJuxtacortLesions numJuxtacortLesions];
    stats.PctJuxtacortLesions = [stats.PctJuxtacortLesions pctJuxtacortLesions];
    
    % Juxtacortical lesions are done by percentage of number of lesions that juxtacortical
    if pctJuxtacortLesions>=juxtaPctThreshRange(2)  % percent of lesions within that distance
        juxtacortical = [juxtacortical ; "yes"];
    elseif pctJuxtacortLesions<=juxtaPctThreshRange(1)
        juxtacortical = [juxtacortical ; "no"];
    else
        juxtacortical = [juxtacortical ; "NA"];  % unsure if juxtacortical because in-between value
    end
    
    % Symmetry
    curSym = summaryInput.LRSymmetry(summaryInput.ID==curACC);
    if abs(curSym)>=symmetryThresh
        symmetry = [symmetry ; "asymmetric"];
    else
        symmetry = [symmetry ; "symmetric"];
    end
    
    realInds=indivLesions.Size_of_Lesion_mm3>20;  % do a cluster threshold here, don't believe tiny lesions
    realLesions=indivLesions.Size_of_Lesion_mm3(realInds);
    clear realInds
    curSize85 = prctile(realLesions,85);
    curMaxSize = max(indivLesions.Size_of_Lesion_mm3);
    stats.MaxLesionSize = [stats.MaxLesionSize max(indivLesions.Size_of_Lesion_mm3)];
    stats.MinLesionSize = [stats.MinLesionSize min(indivLesions.Size_of_Lesion_mm3)];
    stats.MeanLesionSize = [stats.MeanLesionSize mean(indivLesions.Size_of_Lesion_mm3)];
    stats.MedianLesionSize = [stats.MedianLesionSize median(indivLesions.Size_of_Lesion_mm3)];
    stats.LesionSize85pct = [stats.LesionSize85pct curSize85];
    
    % Lesion size empirically done by largest lesion, but keep track of
    % mean, median, min, max, and 85th percentile for stats purposes
    if curMaxSize>=sizeThresh(2)
        lesionSize = [lesionSize ; "large"];
    elseif curMaxSize<=sizeThresh(1)
        lesionSize = [lesionSize ; "small"];
    else
        lesionSize = [lesionSize ; "medium"];
    end
    
    allSizes = indivLesions.Size_of_Lesion_mm3;
    realSizes = allSizes(allSizes>clusterThresh); % lesions must be greater than this number (in mm^3) to count
    clear allSizes
    clear tmpThresh
    curNumLesions = size(realSizes,1);
    stats.NumLesions = [stats.NumLesions curNumLesions];
    
    % number of lesions
    if curNumLesions >= numberThreshRange(2)
        number = [number ; "multiple"];  % should have been named "many" rather than "multiple"
    elseif curNumLesions <= numberThreshRange(1) && curNumLesions >=0
        number = [number ; "single"];  % should have been named "few" rather than "single"
    else
        number = [number ; "NA"];  % in between values, could be noise or could be real, so won't use this feature because unsure
    end
    
    % anterior temporal lobe
    curAntTempVox_L = summaryInput.L_hem_AntTempLobe_vox(summaryInput.ID==curACC);
    curAntTempVox_R = summaryInput.R_hem_AntTempLobe_vox(summaryInput.ID==curACC);
    curAntTempVoxTotal = curAntTempVox_L+curAntTempVox_R;
    curAntTempSym = (curAntTempVox_L-curAntTempVox_R)/(curAntTempVox_L+curAntTempVox_R);
    stats.antTempL = [stats.antTempL curAntTempVox_L];
    stats.antTempR = [stats.antTempR curAntTempVox_R];
    stats.antTempTot = [stats.antTempTot curAntTempVoxTotal];
    stats.antTempSym = [stats.antTempSym curAntTempSym];
    if curAntTempVoxTotal>=antTempThresh && abs(curAntTempSym)<antTempAsymmetryThresh
        anterior_temporal = [anterior_temporal; "yes"];
    else
        anterior_temporal = [anterior_temporal; "no"];
    end
    
    % Mass effect
    curMasseffectRatio = summaryInput.MassEffectAsymRatio(summaryInput.ID==curACC);
    curMasseffectNormed = curMasseffectRatio/curNumLesions; % normed by number of lesions (as calculated above)
    stats.masseffect = [stats.masseffect curMasseffectRatio];
    stats.masseffectNumberNorm = [stats.masseffectNumberNorm curMasseffectNormed];
    if curMasseffectNormed>=masseffectNumberNormThresh(2)
        mass_effect = [mass_effect ; "yes"];
    elseif curMasseffectNormed<=masseffectNumberNormThresh(1)
        mass_effect = [mass_effect ; "no"];
    else  % not sure with intermediate values
        mass_effect = [mass_effect ; "NA"];
    end
    
    % Lobar involvement
    curLobarPercents = [summaryInput.Percentage_Volume_Frontal_Lobe(summaryInput.ID==curACC) ...
        summaryInput.Percentage_Volume_Volume_Parietal_Lobe(summaryInput.ID==curACC) ...
        summaryInput.Percentage_Volume_Volume_Temporal_Lobe(summaryInput.ID==curACC) ...
        summaryInput.Percentage_Volume_Volume_Occipital_Lobe(summaryInput.ID==curACC)];
    stats.LobarPercents = [stats.LobarPercents; curLobarPercents];
    [curLobesSorted, curLobesInds] = sort(curLobarPercents,'descend');
    lobeLabels = ["frontal" "parietal" "temporal" "occipital"]; %indices correspond to these labels, in order
    
    if curLobarPercents(1)>=frontalThresh && curLobarPercents(2)<parietalThresh && curLobarPercents(3)<temporalThresh && curLobarPercents(4)<occipitalThresh % if in 2 lobes, don't count, must be predominantly 1 lobe
        lobar_distribution = [lobar_distribution; "frontal"];
    elseif curLobarPercents(1)<frontalThresh && curLobarPercents(2)>=parietalThresh && curLobarPercents(3)<temporalThresh && curLobarPercents(4)<occipitalThresh % if in 2 lobes, don't count, must be predominantly 1 lobe
        lobar_distribution = [lobar_distribution; "parietal"];
    elseif curLobarPercents(1)<frontalThresh && curLobarPercents(2)<parietalThresh && curLobarPercents(3)>=temporalThresh && curLobarPercents(4)<occipitalThresh % if in 2 lobes, don't count, must be predominantly 1 lobe
        lobar_distribution = [lobar_distribution; "temporal"];
    elseif curLobarPercents(1)<frontalThresh && curLobarPercents(2)<parietalThresh && curLobarPercents(3)<temporalThresh && curLobarPercents(4)>=occipitalThresh % if in 2 lobes, don't count, must be predominantly 1 lobe
        lobar_distribution = [lobar_distribution; "occipital"];
    else
        lobar_distribution = [lobar_distribution; "NA"];
    end
    
    % Ventricular Volume (enlarged or normal)
    curVentVol = summaryInput.Ventricular_Volume_mm3(summaryInput.ID==curACC);
    stats.VentVol = [stats.VentVol; curVentVol];
    
    if curVentVol >= ventThreshRange(2)
        ventVol = [ventVol; "enlarged"];
    elseif curVentVol <= ventThreshRange(1)
        ventVol = [ventVol; "normal"];
    else
        ventVol = [ventVol; "NA"];  % in between volumes, unsure
    end
    
    
    % Lesion Extent (extensive or limited)
    curLesionVol = summaryInput.Total_Size_of_Lesions_mm3(summaryInput.ID==curACC);
    curWMvol = summaryInput.Total_WM_vol_mm3(summaryInput.ID==curACC);
    curLesionExtent = curLesionVol;
    stats.TotLesionVol = [stats.TotLesionVol; curLesionVol];
    stats.WMvol = [stats.WMvol; curWMvol];
    
    
    if curLesionExtent >= extentThreshRange(2)
        lesionExtent = [lesionExtent; "extensive"];
    elseif curLesionExtent <= extentThreshRange(1)
        lesionExtent = [lesionExtent; "limited"];
    else
        lesionExtent = [lesionExtent; "NA"];  % in between values, unsure
    end
end

%% WRITE OUT A CSV FILE WITH ALL THE APPROPRIATE CATEGORIES (BAYES TABLE INPUT FILE)

outputTable = table(acc,age,gender,chronicity,immunocompromised,prodrome,...
    diffusion,enhancement,FLAIR,susceptibility,T1,T2,...
    corpusCallosum,cortex,lobar_distribution,mass_effect,number,periventricular,juxtacortical,...
    lesionSize,symmetry,anterior_temporal,ventVol,lesionExtent,enhancementRatio);


% table needs the output name 'size', but don't want to interfere with the matlab function
outputTable.Properties.VariableNames{'lesionSize'} = 'size';

if saveOutputs
    writetable(outputTable,BayesFile);
    
    save('stats.mat','-struct','stats');
    
    % to load back in, do stats = load('stats.mat');
end

end



