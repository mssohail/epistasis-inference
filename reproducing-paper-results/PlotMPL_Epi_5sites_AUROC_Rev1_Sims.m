%% % code for 5-site simulation Figure 4 run after analysis code and before pltting Figure 4 

clc
% close all
clear all
init_ColorNames_Settings;
% set(0,'DefaultAxesFontName','CMU Serif Roman')
% set(0,'DefaultTextFontName','CMU Serif Roman')
% set(0,'DefaultAxesFontSize',15)
% set(0,'DefaultTextFontSize',15)
warning off
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',7)
set(0,'DefaultTextFontSize',7)


% 0: accessibility defined by rref, i.e., linearly independent col of Hmat
% 1: accessibility defined by polymorphism of site or pair frequency



exportData2Excel = 0;
saveFigs = 0;
% Lin = 10;
thisSet = 1068243%[1068141 1068144 1068241 1068244 1068143 1068142 1068243 1068242]
fileNameContainingDirPath = 'dirNames.txt';
postFiltering = 1; % 0 : pre filtering (out put of AnalysisMPL_filt)
                   % 1 : post filtering
similarityThresh = 0.1; % threshold that sets allowed difference between two numbers to be still declared as same value, ie., 2 = 2.01
noiseThresh = 0.0;
getSysParam;
Tstart = 1;
Tused = 100;%300;
ng = 100;
dT = 10;
numStrainsInInitialPop = 5%20%10;%5; % 5, 10, 30
useStrongReg = 0;

dataFilesWith2Regs = 1;
if(dataFilesWith2Regs == 1)
    regStr1 = '1';
    regStr2 = '1';
end

numItr = 1000;%150%1000%90%90%250;

corrFreqThreshAll = [0 0.01 0.02 0.025 0.03 0.04 0.05 0.075 0.1]; %corrFreqThreshAll = [0 1 2 3 4 5 10];

% Tend = Tused + Tstart - 1;
Tend = Tused + Tstart;

actualT = T/1000;
posSelc = selVal/2/Nin;
negSelc = delSelVal/2/Nin;
lineCol = {'r.-', 'b.-', 'k.-', 'g.-'};
numTps = Tused/dT;
numParam = Lin*(Lin+1)/2;
step = 0.002;
edges = [-10 -0.3:step:0.3 10] + step/2;

maskTri = zeros(Lin, Lin);
for l = 1:Lin
    maskTri(l,l:Lin) = ones(1,Lin-l+1);
end

dirNameScriptFile = pwd;    
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    disp('Error: system si not unix and not PC...')
    pause
end

if(classesOfSites == 2)
    selTypeName = 'PosSel';
elseif(classesOfSites == 3)
    selTypeName = 'PosDelSel';
end
[~, dirNameAnalysis, dirNameFigures] = loadDirNames(fileNameContainingDirPath);
dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
dirNameFigures = [dirNameFigures 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
if(exist(dirNameFigures, 'dir') == 0)
    mkdir(dirNameFigures)        
end
selcCoeffNameColumnOneItr = ['s1  ';'s2  ';'s3  ';'s4  ';'s5  ';'ss12';'ss13';'ss14';'ss15';'ss23';'ss24';'ss25';'ss34';'ss35';'ss45'];
methodCol = [repmat('MPLE     ', 15,1) ; repmat('MPLE-NoMu', 15,1) ; repmat('MPL      ', 15,1) ; repmat('SL       ', 15,1)];

selcCoeffNameAllItr = repmat(' ', numItr*60, 4);
estimatesAllItr = repmat(' ', numItr*60, numParam);
methodsAllItr = repmat(' ', numItr*60, 9);
validAllItr = repmat(' ', numItr*60, 1);

trueClassROC_AllItr = [];
trueClassROCSelc_AllItr = [];
trueClassROCEpi_AllItr = [];
trueClassROCEpiSelc_AllItr = [];

sigmaEstOutLinkEpi_AllItr = [];
sigmaEstOutLinkEpiSelc_AllItr = [];
sigmaEstOutLinkNoMuEpi_AllItr = [];
sigmaEstOutLinkNoMuEpiSelc_AllItr = [];
sigmaEstOutLinkEpiNoE_AllItr = [];
sigmaEstOutLinkEpiNoESelc_AllItr = [];
sigmaEstOutLink_AllItr = [];
sigmaEstOutLinkSelc_AllItr = [];
sigmaEstOutUnLink_AllItr = [];
sigmaEstOutUnLinkSelc_AllItr = [];
perSiteSelction_AllItr  = [];
perSiteSelctionAllEpiTerms_AllItr  = [];

siteToExculeCuzOfLimitedPolyTps_AllItr = [];

sigmaEstOutLinkItr = -1*ones(numItr, 5);
sigmaEstOutLinkEpiItr = -1*ones(numItr, 15);
sigmaEstOutLinkEpiWithRItr = -1*ones(numItr, 15);

aucLinkItr = -1*ones(numItr,2);
aucUnLinkItr = -1*ones(numItr,2);
aucLinkEpiItr = -1*ones(numItr,2);
aucLinkNoMuEpiItr = -1*ones(numItr,2);
aucLinkEpiNoEItr = -1*ones(numItr,2);
aucLinkEpiOffTermsItr = -1*ones(numItr,2);
aucLinkNoMuEpiOffTermsItr = -1*ones(numItr,2);

aucLinkSelcItr = -1*ones(numItr,2);
aucUnLinkSelcItr = -1*ones(numItr,2);
aucLinkEpiSelcItr = -1*ones(numItr,2);
aucLinkNoMuEpiSelcItr = -1*ones(numItr,2);
aucLinkEpiNoESelcItr = -1*ones(numItr,2);
aucLinkEpiOffTermsSelcItr = -1*ones(numItr,2);
aucLinkNoMuEpiOffTermsSelcItr = -1*ones(numItr,2);

aucLinkEpi_only_si_Itr = -1*ones(numItr,2);
aucLinkEpiWithR_only_si_Itr = -1*ones(numItr,2);
aucLinkEpiWithR_only_sij_Itr = -1*ones(numItr,2);
aucLinkEpi_only_sij_Itr = -1*ones(numItr,2);
aucLinkEpi_only_si_SelcItr = -1*ones(numItr,2);
aucLinkEpi_only_sij_SelcItr = -1*ones(numItr,2);
aucLinkEpiWithR_only_si_SelcItr = -1*ones(numItr,2);
aucLinkEpiWithR_only_sij_SelcItr = -1*ones(numItr,2);

absErr_LinkEpi_itr = -1*ones(numItr,1);
absErr_LinkEpi_only_si_itr = -1*ones(numItr,1);
absErr_LinkEpi_only_sij_itr = -1*ones(numItr,1);
absErr_Link_itr = -1*ones(numItr,1);
absErr_UnLink_itr = -1*ones(numItr,1);

absErr_LinkEpi_Selc_itr = -1*ones(numItr,1);
absErr_LinkEpi_only_si_Selc_itr = -1*ones(numItr,1);
absErr_LinkEpi_only_sij_Selc_itr = -1*ones(numItr,1);
absErr_LinkSelc_itr = -1*ones(numItr,1);
absErr_UnLinkSelc_itr = -1*ones(numItr,1);
absErr_LinkEpi_cluster_itr = -1*ones(numItr,1);
absErr_LinkEpi_clusterComb_itr = -1*ones(numItr,1);
absErr_LinkEpi_clusterOnlyWell_itr = -1*ones(numItr,1);
absErr_LinkEpi_clusterOnlyAmb_itr = -1*ones(numItr,1);
absErr_LinkEpi_clusterCombOnlyWell_itr = -1*ones(numItr,1);
absErr_LinkEpi_clusterCombOnlyAmb_itr = -1*ones(numItr,1);

numTerms_LinkEpi_cluster_itr  = -1*ones(numItr,1);
sigmaEstOutLinkEpi_AccAsSum_itr = -1*ones(numItr,1);
rhoItr = zeros(numItr,1);
numSelcSitesItr = zeros(1, numItr);
numSelcSitesEpiItr = zeros(1, numItr);
numSelcSitesEpiItr_si = zeros(1, numItr);
numSelcSitesEpiItr_sij = zeros(1, numItr);
numInAccessible_si = zeros(1, numItr);
numInAccessible_sij = zeros(1, numItr);
numAccessibleAsASum_si = zeros(1, numItr);
numAccessibleAsASum_sij = zeros(1, numItr);
    
allEstEpiItr = zeros(numItr, numParam);
allEstNoMuEpiItr = zeros(numItr, numParam);
allEstEpiNoEItr = zeros(numItr, numParam);
allEstItr = zeros(numItr, Lin);

allSelcSites = zeros(numItr, Lin);
allSelcSitesEpi = zeros(numItr, numParam);
allaccessibleAsSumSitesEpi = zeros(numItr, numParam);

rhoAll_Spearman = -100*ones(numItr, length(corrFreqThreshAll));
pValAll_Spearman = -100*ones(numItr, length(corrFreqThreshAll));
rhoAll_Pearson = -100*ones(numItr, length(corrFreqThreshAll));
pValAll_Pearson = -100*ones(numItr, length(corrFreqThreshAll));

rhoAll_Spearman_MPL = -100*ones(numItr, length(corrFreqThreshAll));
pValAll_Spearman_MPL = -100*ones(numItr, length(corrFreqThreshAll));
rhoAll_Pearson_MPL = -100*ones(numItr, length(corrFreqThreshAll));
pValAll_Pearson_MPL = -100*ones(numItr, length(corrFreqThreshAll));

allStrains32 = ([[zeros(16,1); ones(16, 1)] repmat([zeros(8, 1); ones(8, 1)], 2,1) repmat([0; 0; 0; 0; 1; 1; 1; 1], 4,1) repmat([0; 0; 1; 1], 8,1) repmat([0; 1], 16,1)]);
rhoAll_Spearman32 = -100*ones(numItr, length(corrFreqThreshAll));
pValAll_Spearman32 = -100*ones(numItr, length(corrFreqThreshAll));
rhoAll_Pearson32 = -100*ones(numItr, length(corrFreqThreshAll));
pValAll_Pearson32 = -100*ones(numItr, length(corrFreqThreshAll));

NRMSE_LinkEpi_all = zeros(1, numItr);
NRMSE_LinkEpi_selcSitesEpi = zeros(1, numItr);
NRMSE_LinkEpi_NonSelcSitesEpi = zeros(1, numItr);
NRMSE_LinkEpi_accessibleAsSumSitesEpi = zeros(1, numItr);
NRMSE_LinkEpi_accessibleAsSumSitesEpi_sum = zeros(1, numItr);

RMSE_LinkEpi_all = zeros(1, numItr);
RMSE_LinkEpi_selcSitesEpi = zeros(1, numItr);
RMSE_LinkEpi_NonSelcSitesEpi = zeros(1, numItr);
RMSE_LinkEpi_accessibleAsSumSitesEpi = zeros(1, numItr);
RMSE_LinkEpi_accessibleAsSumSitesEpi_sum = zeros(1, numItr);
RMSE_LinkEpi_inaccessibleSitesEpi = zeros(1, numItr);

trueSelection_selcSitesEpi = 0*ones(1, numItr);
trueSelection_NonSelcSitesEpi = 0*ones(1, numItr);
trueSelection_accessibleAsSumSitesEpi = 0*ones(1, numItr);
trueSelection_accessibleAsSumSitesEpi_sum = 0*ones(1, numItr);
trueSelection_inaccessibleSitesEpi = 0*ones(1, numItr);

epiTermsAccessible = zeros(1, numItr);

timeMPLItr = zeros(1, numItr);
absErr_LinkEpi_all = zeros(1, numItr);
absErr_LinkEpi_selcSitesEpi = zeros(1, numItr);
absErr_LinkEpi_NonSelcSitesEpi = zeros(1, numItr);
condHmatWellItr = zeros(1, numItr);
cond1HmatWellItr = zeros(1, numItr);
rhoItrPear = zeros(1, numItr);
numElementsItrPear = zeros(1, numItr);
rhoItrPear2 = zeros(1, numItr);
numElementsItrPear2 = zeros(1, numItr);
rhoItrSpe = zeros(1, numItr);
numElementsItrSpe = zeros(1, numItr);
rhoItrSpe2 = zeros(1, numItr);
numElementsItrSpe2 = zeros(1, numItr);

diffInSC = zeros(numItr, Lin);
for itr = 1:numItr
    itr
    selcSites = [];
    if(recombination == 1)
        if(dataFilesWith2Regs == 1)
            if(Tstart == 1)
                fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
            else        
                fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
            end        
        else
            if(Tstart == 1)
                fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];

                fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
            end
        end
        
        if(useStrongReg == 1)
            fileName = [fileName(1:end-4) '_useStrong.mat'];
        end
        load([dirNameAnalysis fileName], 'L', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
             'sigmaEstOutLink', 'sigmaEstOutLinkEpi', 'sigmaEstOutLinkEpiWithR', 'sigmaEstOutLinkNoMuEpi',...
             'sigmaEstOutUnLink', 'sigmaEstOutLinkEpiNoE', 'sumAijWithLink', ...
             'q', 'q11', 'priorConst', 'timeWholeCodeMPL',...'perSiteSelctionSelc', 'selcSites', ... 'sigmaEstOutLinkSelc', 'sigmaEstOutLinkNoMuSelc','sigmaEstOutUnLinkSelc', 'sigmaEstOutUnLinkNoMuSelc', ...
             'deltaLogLikeli', 'deltaLogLikeli50', 'deltaLogLikeli100',...
             'deltaLogLikeliMPL', 'deltaLogLikeliSL', 'Hmat', ...
             'masterStrainList', 'currentNumOfCirStrainsInAllPop', 'masterSelectionFactor', 'freqStatsAllStrains', ...
             'Reg2Mat', 'qExt_tLast', 'qExt_tStart', 'vExt_EstOutLink');
         
    else
        if(dataFilesWith2Regs == 1)
            if(Tstart == 1)
                fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
                
            else        
                fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
            end        
        else
            if(Tstart == 1)
                fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
            else        
                fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
            end
        end
        load([dirNameAnalysis fileName], 'Lin', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
             'sigmaEstOutLink', 'sigmaEstOutLinkEpi', 'sigmaEstOutLinkNoMuEpi',...
             'sigmaEstOutUnLink', 'sigmaEstOutLinkEpiNoE', 'sumAijWithLink', ...
             'q', 'q11', 'priorConst', 'timeWholeCodeMPL',...'perSiteSelctionSelc', 'selcSites', ... 'sigmaEstOutLinkSelc', 'sigmaEstOutLinkNoMuSelc','sigmaEstOutUnLinkSelc', 'sigmaEstOutUnLinkNoMuSelc', ...
             'deltaLogLikeli', 'deltaLogLikeli50', 'deltaLogLikeli100',...
             'deltaLogLikeliMPL', 'deltaLogLikeliSL', 'Hmat', ...
             'masterStrainList', 'currentNumOfCirStrainsInAllPop', 'masterSelectionFactor', 'freqStatsAllStrains', ...
             'Reg2Mat', 'qExt_tLast', 'qExt_tStart', 'vExt_EstOutLink');
         
    end
%     [selcSitesEpi, wellCondSites, illCondSites, cluster] = findIndColsOfHmat(Hmat)
%     pause
   
    
    trueClassROC = 3*ones(1, Lin);
    trueClassROC(perSiteSelction == posSelc) = 4;
    trueClassROC(perSiteSelction == negSelc) = 5;

    sigmaEstOutLinkItr(itr,:) = sigmaEstOutLink';
    sigmaEstOutLinkEpiItr(itr,:) = sigmaEstOutLinkEpi';
    sigmaEstOutLinkEpiWithRItr(itr,:) = sigmaEstOutLinkEpiWithR';
    
    perSiteSelctionAllEpiTerms = [];
    for l = 1:Lin
        perSiteSelctionAllEpiTerms = [perSiteSelctionAllEpiTerms perSiteSelctionEpi(l,l+1:end)];
    end

    perSiteSelctionAllEpiTerms = [perSiteSelction perSiteSelctionAllEpiTerms];

    trueClassROCEpi = 3*ones(1, (Lin*(Lin + 1)/2));
    trueClassROCEpi(perSiteSelction == posSelc) = 4;
    trueClassROCEpi(perSiteSelction == negSelc) = 5;
    
    [selcSitesEpi1, wellCondSites1, illCondSites1, cluster] = findIndColsOfHmat(Hmat*dT);
    
    
    selcSitesEpi = selcSitesEpi1;
    
    wellCondSites = wellCondSites1;
    illCondSites = setdiff(1:numParam,wellCondSites);
    
    sitesWithZeroFreq = find(sum(Hmat) == 0);
    accessibleAsSumSites = setdiff(illCondSites, sitesWithZeroFreq);
    
%      if(itr <= 10)
%     [selcSitesEpi, wellCondSites, illCondSites, cluster] = findIndColsOfHmat(Hmat)
%    
%      
%      newEst = (Hmat(selcSitesEpi,selcSitesEpi)*dT + Reg2Mat(selcSitesEpi,selcSitesEpi))\(qExt_tLast(selcSitesEpi) - qExt_tStart(selcSitesEpi) - vExt_EstOutLink(selcSitesEpi)')';
%      oldEst = sigmaEstOutLinkEpi(selcSitesEpi);
%      trueEst = perSiteSelctionAllEpiTerms(selcSitesEpi);
%      disp('newEst trueEst oldEst wellCondSitesIndex')
%      [newEst trueEst' oldEst wellCondSites']
%       pause
%      end
    
    
    zeroCol = sum(Hmat == 0) == (Lin*(Lin+1)/2);
    
    epiTermsAccessible(itr) = sum(selcSitesEpi(Lin+1:end));
    if(postFiltering == 1)
        
        selcSites = findIndColsOfHmat(Hmat(1:Lin,1:Lin)*dT);
        
        numSelcSitesItr(itr) = sum(selcSites);
        numSelcSitesEpiItr(itr) = sum(selcSitesEpi);
        numSelcSitesEpiItr_si(itr) = sum(selcSitesEpi(1:Lin));
        numSelcSitesEpiItr_sij(itr) = sum(selcSitesEpi(Lin+1:end));
        numInAccessible_si(itr) = sum(zeroCol(1:Lin));
        numInAccessible_sij(itr) = sum(zeroCol(Lin+1:(Lin*(Lin+1)/2)));
        numAccessibleAsASum_si(itr) = Lin - numInAccessible_si(itr) - numSelcSitesEpiItr_si(itr);
        numAccessibleAsASum_sij(itr) = (Lin*(Lin-1)/2)  - numInAccessible_sij(itr) - numSelcSitesEpiItr_sij(itr);
        
        perSiteSelctionSelc = perSiteSelction(selcSites);
        perSiteSelctionAllEpiTermsSelc = perSiteSelctionAllEpiTerms(selcSitesEpi);
        trueClassROCSelc = trueClassROC(selcSites);
        trueClassROCEpiSelc = trueClassROCEpi(selcSitesEpi);
        
        
        sigmaEstOutLinkEpiSelc = sigmaEstOutLinkEpi(selcSitesEpi);
        sigmaEstOutLinkEpiWithRSelc = sigmaEstOutLinkEpiWithR(selcSitesEpi);
        sigmaEstOutLinkNoMuEpiSelc = sigmaEstOutLinkNoMuEpi(selcSitesEpi);
        sigmaEstOutLinkEpiNoESelc = sigmaEstOutLinkEpiNoE(selcSitesEpi);
        sigmaEstOutLinkSelc = sigmaEstOutLink(selcSites);
        sigmaEstOutUnLinkSelc = sigmaEstOutUnLink(selcSites);        
    end
    allSelcSites(itr,:) = selcSites;
    allSelcSitesEpi(itr,:) = selcSitesEpi;
    
    timeMPLItr(itr) = timeWholeCodeMPL(itr);
    
    diffInSC(itr,:) = abs(sigmaEstOutLinkEpiWithR(1:Lin) - sigmaEstOutLink);
    
    % calculate AUROC
    %--------------------------------------------------
    neutralLimit = 0.005;
    posOnlyItrTemp = perSiteSelction > neutralLimit;%== max(perSiteSelction);
    negOnlyItrTemp = perSiteSelction < -neutralLimit;%== min(perSiteSelction);
    posOnlySelcItrTemp = perSiteSelctionSelc > neutralLimit;%== max(perSiteSelctionSelc);
    negOnlySelcItrTemp = perSiteSelctionSelc < -neutralLimit;%== min(perSiteSelctionSelc);
    [aucLinkItrTemp, aucLinkSelcItrTemp] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutLink, sigmaEstOutLinkSelc, posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp);
    [aucUnLinkItrTemp, aucUnLinkSelcItrTemp] = calcAUROC(perSiteSelction, perSiteSelctionSelc, sigmaEstOutUnLink, sigmaEstOutUnLinkSelc, posOnlyItrTemp, negOnlyItrTemp, posOnlySelcItrTemp, negOnlySelcItrTemp);
    
    % all EPi terms
    posOnlyItrTemp_AllEpiTerms = perSiteSelctionAllEpiTerms > neutralLimit;%== max(perSiteSelction);
    negOnlyItrTemp_AllEpiTerms = perSiteSelctionAllEpiTerms < -neutralLimit;%== min(perSiteSelction);
    posOnlySelcItrTemp_AllEpiTerms = perSiteSelctionAllEpiTermsSelc > neutralLimit;%== max(perSiteSelctionSelc);
    negOnlySelcItrTemp_AllEpiTerms = perSiteSelctionAllEpiTermsSelc < -neutralLimit;%== min(perSiteSelctionSelc);
    [aucLinkEpiItrTemp, aucLinkEpiSelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkEpi, sigmaEstOutLinkEpiSelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);
    [aucLinkEpiWithRItrTemp, aucLinkEpiWithRSelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkEpiWithR, sigmaEstOutLinkEpiWithRSelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);
    [aucLinkNoMuEpiItrTemp, aucLinkNoMuEpiSelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkNoMuEpi, sigmaEstOutLinkNoMuEpiSelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);
    [aucLinkEpiNoEItrTemp, aucLinkEpiNoESelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkEpiNoE, sigmaEstOutLinkEpiNoESelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);

    % only s_i term sof Epi
    sitesEpi_only_si = logical([ones(1, Lin), zeros(1, Lin*(Lin-1)/2)]);
    selcSitesEpi_only_si = selcSitesEpi & sitesEpi_only_si;
    perSiteSelction_only_si = perSiteSelctionAllEpiTerms(sitesEpi_only_si);
    perSiteSelction_only_si_Selc = perSiteSelctionAllEpiTerms(selcSitesEpi_only_si);
    sigmaEstOutLinkEpi_only_si = sigmaEstOutLinkEpi(sitesEpi_only_si);
    sigmaEstOutLinkEpi_only_si_Selc = sigmaEstOutLinkEpi(selcSitesEpi_only_si);
    sigmaEstOutLinkEpiWithR_only_si = sigmaEstOutLinkEpiWithR(sitesEpi_only_si);
    sigmaEstOutLinkEpiWithR_only_si_Selc = sigmaEstOutLinkEpiWithR(selcSitesEpi_only_si);
 
    posOnlyItrTemp_only_si = perSiteSelction_only_si > neutralLimit; %posOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
    negOnlyItrTemp_only_si = perSiteSelction_only_si < -neutralLimit; %negOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
    posOnlySelcItrTemp_only_si = perSiteSelction_only_si_Selc > neutralLimit;
    negOnlySelcItrTemp_only_si = perSiteSelction_only_si_Selc < -neutralLimit;
    [aucLinkEpiItrTemp_only_si, aucLinkEpiSelcItrTemp_only_si] = calcAUROC(perSiteSelction_only_si, perSiteSelction_only_si_Selc, sigmaEstOutLinkEpi_only_si, sigmaEstOutLinkEpi_only_si_Selc, posOnlyItrTemp_only_si, negOnlyItrTemp_only_si, posOnlySelcItrTemp_only_si, negOnlySelcItrTemp_only_si);
    [aucLinkEpiWithRItrTemp_only_si, aucLinkEpiWithRSelcItrTemp_only_si] = calcAUROC(perSiteSelction_only_si, perSiteSelction_only_si_Selc, sigmaEstOutLinkEpiWithR_only_si, sigmaEstOutLinkEpiWithR_only_si_Selc, posOnlyItrTemp_only_si, negOnlyItrTemp_only_si, posOnlySelcItrTemp_only_si, negOnlySelcItrTemp_only_si);
 
    % only s_ij term sof Epi
    sitesEpi_only_sij = logical([zeros(1, Lin), ones(1, Lin*(Lin-1)/2)]);
    selcSitesEpi_only_sij = selcSitesEpi & sitesEpi_only_sij;
    perSiteSelction_only_sij = perSiteSelctionAllEpiTerms(sitesEpi_only_sij);
    perSiteSelction_only_sij_Selc = perSiteSelctionAllEpiTerms(selcSitesEpi_only_sij);
    sigmaEstOutLinkEpi_only_sij = sigmaEstOutLinkEpi(sitesEpi_only_sij);
    sigmaEstOutLinkEpi_only_sij_Selc = sigmaEstOutLinkEpi(selcSitesEpi_only_sij);
    sigmaEstOutLinkEpiWithR_only_sij = sigmaEstOutLinkEpiWithR(sitesEpi_only_sij);
    sigmaEstOutLinkEpiWithR_only_sij_Selc = sigmaEstOutLinkEpiWithR(selcSitesEpi_only_sij);
 
    posOnlyItrTemp_only_sij = perSiteSelction_only_sij > neutralLimit; %posOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
    negOnlyItrTemp_only_sij = perSiteSelction_only_sij < -neutralLimit; %negOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
    posOnlySelcItrTemp_only_sij = perSiteSelction_only_sij_Selc > neutralLimit;
    negOnlySelcItrTemp_only_sij = perSiteSelction_only_sij_Selc < -neutralLimit;
    [aucLinkEpiItrTemp_only_sij, aucLinkEpiSelcItrTemp_only_sij] = calcAUROC(perSiteSelction_only_sij, perSiteSelction_only_sij_Selc, sigmaEstOutLinkEpi_only_sij, sigmaEstOutLinkEpi_only_sij_Selc, posOnlyItrTemp_only_sij, negOnlyItrTemp_only_sij, posOnlySelcItrTemp_only_sij, negOnlySelcItrTemp_only_sij);
    [aucLinkEpiWithRItrTemp_only_sij, aucLinkEpiWithRSelcItrTemp_only_sij] = calcAUROC(perSiteSelction_only_sij, perSiteSelction_only_sij_Selc, sigmaEstOutLinkEpiWithR_only_sij, sigmaEstOutLinkEpiWithR_only_sij_Selc, posOnlyItrTemp_only_sij, negOnlyItrTemp_only_sij, posOnlySelcItrTemp_only_sij, negOnlySelcItrTemp_only_sij);
    
    zeroThresh = 0.001;%0.00025;
    naeReg = 0.001;%0.00025;
    
    [absErr_LinkEpi_only_si_Selc_temp, absErr_LinkEpi_only_si_temp, ...
    absErr_LinkEpi_only_sij_Selc_temp, absErr_LinkEpi_only_sij_temp, ...
    absErr_LinkEpi_Selc_temp, absErr_LinkEpi_temp, absErr_LinkSelc_temp, absErr_UnLinkSelc_temp, ...
    absErr_Link_temp, absErr_UnLink_temp, absErr_LinkEpi_cluster_temp, absErr_LinkEpi_clusterOnlyWell_temp, ...
    absErr_LinkEpi_clusterOnlyAmb_temp, absErr_LinkEpi_clusterComb_temp, ...
    absErr_LinkEpi_clusterCombOnlyWell_temp, absErr_LinkEpi_clusterCombOnlyAmb_temp] = calcNAE(itr, zeroThresh, naeReg, ...
            selcSitesEpi_only_si, perSiteSelction_only_si_Selc, sigmaEstOutLinkEpi_only_si_Selc, ...
            perSiteSelction_only_si, sigmaEstOutLinkEpi_only_si,  selcSitesEpi_only_sij, perSiteSelction_only_sij_Selc, ...
            sigmaEstOutLinkEpi_only_sij_Selc, perSiteSelction_only_sij, sigmaEstOutLinkEpi_only_sij, selcSitesEpi, ...
            perSiteSelctionAllEpiTerms, sigmaEstOutLinkEpiSelc, sigmaEstOutLinkEpi, selcSites, perSiteSelction, ...
            sigmaEstOutLinkSelc, sigmaEstOutUnLinkSelc, sigmaEstOutLink, sigmaEstOutUnLink, wellCondSites, cluster);
    
    absErr_LinkEpi_all(itr) = sum(abs(perSiteSelctionAllEpiTerms - sigmaEstOutLinkEpi'))/(sum(abs(perSiteSelctionAllEpiTerms)));
    NRMSE_LinkEpi_all(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms - sigmaEstOutLinkEpi').^2)/sum(abs(perSiteSelctionAllEpiTerms).^2));
    RMSE_LinkEpi_all(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms - sigmaEstOutLinkEpi').^2))/15;
    if(sum(selcSitesEpi) == 0)
        absErr_LinkEpi_selcSitesEpi(itr) = -1;
        NRMSE_LinkEpi_selcSitesEpi(itr) = -1;
        RMSE_LinkEpi_selcSitesEpi(itr) = -1;
    else
        absErr_LinkEpi_selcSitesEpi(itr) = sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi) - sigmaEstOutLinkEpi(selcSitesEpi)'))/(sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi))));
        NRMSE_LinkEpi_selcSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi) - sigmaEstOutLinkEpi(selcSitesEpi)').^2)/sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi)).^2));
        RMSE_LinkEpi_selcSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi) - sigmaEstOutLinkEpi(selcSitesEpi)').^2))/sum(selcSitesEpi);
    end
    
    if(sum(selcSitesEpi) == Lin*(Lin+1)/2)
        absErr_LinkEpi_NonSelcSitesEpi(itr) = -1;
        NRMSE_LinkEpi_NonSelcSitesEpi(itr) = -1;
        NRMSE_LinkEpi_accessibleAsSumSitesEpi(itr) = -1;
        NRMSE_LinkEpi_accessibleAsSumSitesEpi_sum(itr) = -1;
        RMSE_LinkEpi_NonSelcSitesEpi(itr) = -1;
        sigmaEstOutLinkEpi_AccAsSum_itr(itr) = 0;
        RMSE_LinkEpi_accessibleAsSumSitesEpi(itr) = -1;
        RMSE_LinkEpi_accessibleAsSumSitesEpi_sum(itr) = -1;
    else
        absErr_LinkEpi_NonSelcSitesEpi(itr) = sum(abs(perSiteSelctionAllEpiTerms(~selcSitesEpi) - sigmaEstOutLinkEpi(~selcSitesEpi)'))/(sum(abs(perSiteSelctionAllEpiTerms(~selcSitesEpi))));
        NRMSE_LinkEpi_NonSelcSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(~selcSitesEpi) - sigmaEstOutLinkEpi(~selcSitesEpi)').^2)/sum(abs(perSiteSelctionAllEpiTerms(~selcSitesEpi)).^2));
        NRMSE_LinkEpi_accessibleAsSumSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(accessibleAsSumSites) - sigmaEstOutLinkEpi(accessibleAsSumSites)').^2)/sum(abs(perSiteSelctionAllEpiTerms(accessibleAsSumSites)).^2));
        NRMSE_LinkEpi_accessibleAsSumSitesEpi_sum(itr) = sqrt((abs(sum(perSiteSelctionAllEpiTerms(accessibleAsSumSites)) - sum(sigmaEstOutLinkEpi(accessibleAsSumSites))').^2)/(abs(sum(perSiteSelctionAllEpiTerms(accessibleAsSumSites))).^2));
        sigmaEstOutLinkEpi_AccAsSum_itr(itr) = sum(sigmaEstOutLinkEpi(accessibleAsSumSites));
        RMSE_LinkEpi_NonSelcSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(~selcSitesEpi) - sigmaEstOutLinkEpi(~selcSitesEpi)').^2))/sum(~selcSitesEpi);
        RMSE_LinkEpi_accessibleAsSumSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(accessibleAsSumSites) - sigmaEstOutLinkEpi(accessibleAsSumSites)').^2))/length(accessibleAsSumSites);
        RMSE_LinkEpi_accessibleAsSumSitesEpi_sum(itr) = sqrt((abs(sum(perSiteSelctionAllEpiTerms(accessibleAsSumSites)) - sum(sigmaEstOutLinkEpi(accessibleAsSumSites))').^2))/length(accessibleAsSumSites);
    end
    
    
    if(isempty(sitesWithZeroFreq))
        RMSE_LinkEpi_inaccessibleSitesEpi(itr) = -1;
    else
        RMSE_LinkEpi_inaccessibleSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(sitesWithZeroFreq) - sigmaEstOutLinkEpi(sitesWithZeroFreq)').^2))/length(sitesWithZeroFreq);
    end
    temp1 = zeros(1,15);
    temp1(accessibleAsSumSites) = 1;
    allaccessibleAsSumSitesEpi(itr,:) = temp1;
    
    trueSelection_selcSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi)).^2));%sum(perSiteSelctionAllEpiTerms(selcSitesEpi));
    trueSelection_NonSelcSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(~selcSitesEpi)).^2));%sum(perSiteSelctionAllEpiTerms(~selcSitesEpi));
    trueSelection_accessibleAsSumSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(accessibleAsSumSites)).^2));%sum(perSiteSelctionAllEpiTerms(accessibleAsSumSites));
    trueSelection_accessibleAsSumSitesEpi_sum(itr) = sqrt((abs(sum(perSiteSelctionAllEpiTerms(accessibleAsSumSites))).^2));%sum(perSiteSelctionAllEpiTerms(accessibleAsSumSites));
    trueSelection_inaccessibleSitesEpi(itr) = sqrt(sum(abs(perSiteSelctionAllEpiTerms(sitesWithZeroFreq)).^2));%sum(perSiteSelctionAllEpiTerms(accessibleAsSumSites));
    
       
    absErr_LinkEpi_only_si_Selc_itr(itr) = absErr_LinkEpi_only_si_Selc_temp;
    absErr_LinkEpi_only_si_itr(itr) = absErr_LinkEpi_only_si_temp;
    absErr_LinkEpi_only_sij_Selc_itr(itr) = absErr_LinkEpi_only_sij_Selc_temp;
    absErr_LinkEpi_only_sij_itr(itr) = absErr_LinkEpi_only_sij_temp;
    absErr_LinkEpi_Selc_itr(itr) = absErr_LinkEpi_Selc_temp;
    absErr_LinkEpi_itr(itr) = absErr_LinkEpi_temp;
    absErr_LinkSelc_itr(itr) = absErr_LinkSelc_temp;
    absErr_UnLinkSelc_itr(itr) = absErr_UnLinkSelc_temp;
    absErr_Link_itr(itr) = absErr_Link_temp;
    absErr_UnLink_itr(itr) = absErr_UnLink_temp;
    absErr_LinkEpi_cluster_itr(itr) = absErr_LinkEpi_cluster_temp;
    absErr_LinkEpi_clusterOnlyWell_itr(itr) = absErr_LinkEpi_clusterOnlyWell_temp;
    absErr_LinkEpi_clusterOnlyAmb_itr(itr) = absErr_LinkEpi_clusterOnlyAmb_temp;
    absErr_LinkEpi_clusterComb_itr(itr) = absErr_LinkEpi_clusterComb_temp;
    absErr_LinkEpi_clusterCombOnlyWell_itr(itr) = absErr_LinkEpi_clusterCombOnlyWell_temp;
    absErr_LinkEpi_clusterCombOnlyAmb_itr(itr) = absErr_LinkEpi_clusterCombOnlyAmb_temp;

    if(absErr_LinkEpi_clusterComb_itr(itr) > 40)
%         pause
    end
    
    % Save all variables for this itr
    %--------------------------------------

    % Save AUROC unfiltered
    aucLinkEpiItr(itr,:) = aucLinkEpiItrTemp;
    aucLinkNoMuEpiItr(itr,:) = aucLinkNoMuEpiItrTemp;
    aucLinkEpiNoEItr(itr,:) = aucLinkEpiNoEItrTemp;
    aucLinkItr(itr,:) = aucLinkItrTemp;
    aucUnLinkItr(itr,:) = aucUnLinkItrTemp;

    % Save AUROC filtered
    aucLinkEpiSelcItr(itr,:) = aucLinkEpiSelcItrTemp;
    aucLinkNoMuEpiSelcItr(itr,:) = aucLinkNoMuEpiSelcItrTemp;
    aucLinkEpiNoESelcItr(itr,:) = aucLinkEpiNoESelcItrTemp;
    aucLinkSelcItr(itr,:) = aucLinkSelcItrTemp;
    aucUnLinkSelcItr(itr,:) = aucUnLinkSelcItrTemp;

    aucLinkEpi_only_si_Itr(itr,:) = aucLinkEpiItrTemp_only_si;
    aucLinkEpi_only_sij_Itr(itr,:) = aucLinkEpiItrTemp_only_sij;
    aucLinkEpiWithR_only_si_Itr(itr,:) = aucLinkEpiWithRItrTemp_only_si;
    aucLinkEpiWithR_only_sij_Itr(itr,:) = aucLinkEpiWithRItrTemp_only_sij;

    aucLinkEpi_only_si_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_si;
    aucLinkEpi_only_sij_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_sij;
    aucLinkEpiWithR_only_si_SelcItr(itr,:) = aucLinkEpiWithRSelcItrTemp_only_si;
    aucLinkEpiWithR_only_sij_SelcItr(itr,:) = aucLinkEpiWithRSelcItrTemp_only_sij;
    
    %[selcSitesEpi, wellCondSites, illCondSites, cluster] = findIndColsOfHmat(Hmat*dT);
    
    condHmatWellItr(itr) = cond(Hmat(wellCondSites,wellCondSites)*dT);
    cond1HmatWellItr(itr) = cond(Hmat(wellCondSites1,wellCondSites1)*dT);
    %cond2HmatWellItr(itr) = cond(Hmat(wellCondSites2,wellCondSites2)*dT);
    
    [rhoItrPear(itr), pval, perSiteSelcTemp80SortedPear, sigmaEstOutLinkEpiTemp80SortedPear] = getCorr(perSiteSelctionEpi, perSiteSelctionAllEpiTerms, cluster, sigmaEstOutLinkEpi, wellCondSites, 'Pearson');
    numElementsItrPear(itr) = length(perSiteSelcTemp80SortedPear);
    [rhoItrPear2(itr), pval, perSiteSelcTemp80SortedPear2, sigmaEstOutLinkEpiTemp80SortedPear2] = getCorr2(perSiteSelctionEpi, perSiteSelctionAllEpiTerms, cluster, sigmaEstOutLinkEpi, wellCondSites, 'Pearson');
    numElementsItrPear2(itr) = length(perSiteSelcTemp80SortedPear2);
    [rhoItrSpe(itr), pval, perSiteSelcTemp80SortedSpe, sigmaEstOutLinkEpiTemp80SortedSpe] = getCorr(perSiteSelctionEpi, perSiteSelctionAllEpiTerms, cluster, sigmaEstOutLinkEpi, wellCondSites, 'Spearman');
    numElementsItrSpe(itr) = length(perSiteSelcTemp80SortedSpe);
    [rhoItrSpe2(itr), pval, perSiteSelcTemp80SortedSpe2, sigmaEstOutLinkEpiTemp80SortedSpe2] = getCorr2(perSiteSelctionEpi, perSiteSelctionAllEpiTerms, cluster, sigmaEstOutLinkEpi, wellCondSites, 'Spearman');
    numElementsItrSpe2(itr) = length(perSiteSelcTemp80SortedSpe2);
    
    % Save estimated SC
    allEstEpiItr(itr,:) = sigmaEstOutLinkEpi';
    allEstNoMuEpiItr(itr,:) = sigmaEstOutLinkNoMuEpi';
    allEstEpiNoEItr(itr,:) = sigmaEstOutLinkEpiNoE';
    allEstItr(itr,:) = sigmaEstOutLink';    
    
    % prep data from export to R
    % --------------------------
    selcCoeffNameAllItr((itr-1)*60 + 1:itr*60,:) = repmat(selcCoeffNameColumnOneItr, 4, 1);
    
    temp15 = num2str([sigmaEstOutLinkEpi; sigmaEstOutLinkNoMuEpi]);
    
    temp20 = num2str(sigmaEstOutLink); 
    temp21 = repmat(' ', 10, size(temp20, 2));
    temp25 = [temp20; temp21];
    
    temp30 = num2str(sigmaEstOutUnLink); 
    temp31 = repmat(' ', 10, size(temp30, 2));
    temp35 = [temp30; temp31];
    
    %[(itr-1)*60+1  itr*60-30  (itr-1)*60+1+30  itr*60-15 (itr-1)*60+1+45   itr*60
    estimatesAllItr((itr-1)*60+1:itr*60-30,1:size(temp15,2)) = temp15;
    estimatesAllItr((itr-1)*60+1+30:itr*60-15,1:size(temp25,2)) = temp25;
    estimatesAllItr((itr-1)*60+1+45:itr*60,1:size(temp35,2)) = temp35;

    methodsAllItr((itr-1)*60 + 1:itr*60,:) = methodCol;
    
    validAllItr((itr-1)*60+1:itr*60-30,1) = num2str([selcSitesEpi selcSitesEpi]');
    validAllItr((itr-1)*60+1+30:itr*60-15,1) = [num2str(selcSites'); repmat(' ', 10, 1)];
    validAllItr((itr-1)*60+1+45:itr*60,1) = [num2str(selcSites'); repmat(' ', 10, 1)];
    
    
    % for scatter plot of haplotype SC
    %------------------------------------
    
    % make upper tri matrix of estimated FL
    for jk = 1:Lin
        if(jk == 1)
            stInd = 1;
            stpInd = stInd + Lin - 1;
            sigmaLinkEpiMtx = diag(sigmaEstOutLinkEpi(stInd:stpInd));
        else
            stInd = stpInd + 1;
            stpInd = stInd + (Lin - jk);
            %[stInd stpInd jk-1 jk]
            sigmaLinkEpiMtx(jk-1,jk:end) = sigmaEstOutLinkEpi(stInd:stpInd);
        end
    end
    
    allStrains = masterStrainList(1:currentNumOfCirStrainsInAllPop,:);
    allStrainsSelFac = masterSelectionFactor(1:currentNumOfCirStrainsInAllPop);
    diffFreq = freqStatsAllStrains(1:currentNumOfCirStrainsInAllPop,2) - freqStatsAllStrains(1:currentNumOfCirStrainsInAllPop,1);
    
    estSelFac = zeros(currentNumOfCirStrainsInAllPop, 1);
    estSelFac_MPL = zeros(currentNumOfCirStrainsInAllPop, 1);
    for jk = 1:currentNumOfCirStrainsInAllPop
        estSelFac(jk) = sum(sum(((allStrains(jk,:)'*allStrains(jk,:))).*maskTri.*sigmaLinkEpiMtx));
        estSelFac_MPL(jk) = allStrains(jk,:)*sigmaEstOutLink;
    end
    
    estSelFac32 = zeros(32, 1);
    allStrainsSelFac32 = zeros(32, 1);
    for jk = 1:32
        estSelFac32(jk) = sum(sum(((allStrains32(jk,:)'*allStrains32(jk,:))).*maskTri.*sigmaLinkEpiMtx));
        allStrainsSelFac32(jk) = sum(sum(((allStrains32(jk,:)'*allStrains32(jk,:))).*maskTri.*perSiteSelctionEpi));
    end
    % calc correlation at threshold frequency of 0, 1, 2, 3, 4, 5, 10
    
    for jk = 1:length(corrFreqThreshAll)
        corrFreqThresh = corrFreqThreshAll(jk);
        qa = diffFreq >  corrFreqThresh;
        if(sum(qa) >=3)
            [rhoTemp,pValTemp] = corr(allStrainsSelFac(qa), estSelFac(qa), 'type', 'Spearman');
            rhoAll_Spearman(itr, jk) = rhoTemp;
            pValAll_Spearman(itr, jk) = pValTemp;
            [rhoTemp,pValTemp] = corr(allStrainsSelFac(qa), estSelFac(qa), 'type', 'Pearson');
            rhoAll_Pearson(itr, jk) = rhoTemp;
            pValAll_Pearson(itr, jk) = pValTemp;
            
            [rhoTemp,pValTemp] = corr(allStrainsSelFac(qa), estSelFac_MPL(qa), 'type', 'Spearman');
            rhoAll_Spearman_MPL(itr, jk) = rhoTemp;
            pValAll_Spearman_MPL(itr, jk) = pValTemp;
            [rhoTemp,pValTemp] = corr(allStrainsSelFac(qa), estSelFac_MPL(qa), 'type', 'Pearson');
            rhoAll_Pearson_MPL(itr, jk) = rhoTemp;
            pValAll_Pearson_MPL(itr, jk) = pValTemp;
        else
        end
    end
    
    for jk = 1:length(corrFreqThreshAll)
        corrFreqThresh = corrFreqThreshAll(jk);
        [rhoTemp,pValTemp] = corr(allStrainsSelFac32, estSelFac32, 'type', 'Spearman');
        rhoAll_Spearman32(itr, jk) = rhoTemp;
        pValAll_Spearman32(itr, jk) = pValTemp;
        [rhoTemp,pValTemp] = corr(allStrainsSelFac32, estSelFac32, 'type', 'Pearson');
        rhoAll_Pearson32(itr, jk) = rhoTemp;
        pValAll_Pearson32(itr, jk) = pValTemp;
    end
end

%%

meanAllNormAbsErrs_Selc = [mean(absErr_LinkEpi_clusterComb_itr(absErr_LinkEpi_clusterComb_itr ~= -1)) ...
           mean(absErr_LinkEpi_cluster_itr(absErr_LinkEpi_cluster_itr ~= -1)) ...
           mean(absErr_LinkEpi_Selc_itr(absErr_LinkEpi_Selc_itr ~= -1)) ...
           mean(absErr_LinkEpi_only_si_Selc_itr(absErr_LinkEpi_only_si_Selc_itr ~= -1)) ...
           mean(absErr_LinkEpi_only_sij_Selc_itr(absErr_LinkEpi_only_sij_Selc_itr ~= -1)) ...
           mean(absErr_LinkSelc_itr(absErr_LinkSelc_itr ~= -1)) ...
           mean(absErr_UnLinkSelc_itr(absErr_UnLinkSelc_itr ~= -1))];

meanAllNormAbsErrs = [mean(absErr_LinkEpi_clusterComb_itr(absErr_LinkEpi_clusterComb_itr ~= -1)) ...
           mean(absErr_LinkEpi_clusterCombOnlyWell_itr(absErr_LinkEpi_clusterCombOnlyWell_itr ~= -1)) ...
           mean(absErr_LinkEpi_clusterCombOnlyAmb_itr(absErr_LinkEpi_clusterCombOnlyAmb_itr ~= -1)) ...
           mean(absErr_LinkEpi_cluster_itr(absErr_LinkEpi_cluster_itr ~= -1)) ...
           mean(absErr_LinkEpi_clusterOnlyWell_itr(absErr_LinkEpi_clusterOnlyWell_itr ~= -1)) ...
           mean(absErr_LinkEpi_clusterOnlyAmb_itr(absErr_LinkEpi_clusterOnlyAmb_itr ~= -1)) ...
           mean(absErr_LinkEpi_itr(absErr_LinkEpi_itr ~= -1)) ...
           mean(absErr_LinkEpi_only_si_itr(absErr_LinkEpi_only_si_itr ~= -1)) ...
           mean(absErr_LinkEpi_only_sij_itr(absErr_LinkEpi_only_sij_itr ~= -1)) ...
           mean(absErr_Link_itr(absErr_Link_itr ~= -1)) ...
           mean(absErr_UnLink_itr(absErr_UnLink_itr ~= -1))];

medianAllNormAbsErrs_Selc = [median(absErr_LinkEpi_clusterComb_itr(absErr_LinkEpi_clusterComb_itr ~= -1)) ...
           median(absErr_LinkEpi_cluster_itr(absErr_LinkEpi_cluster_itr ~= -1)) ...
           median(absErr_LinkEpi_Selc_itr(absErr_LinkEpi_Selc_itr ~= -1)) ...
           median(absErr_LinkEpi_only_si_Selc_itr(absErr_LinkEpi_only_si_Selc_itr ~= -1)) ...
           median(absErr_LinkEpi_only_sij_Selc_itr(absErr_LinkEpi_only_sij_Selc_itr ~= -1)) ...
           median(absErr_LinkSelc_itr(absErr_LinkSelc_itr ~= -1)) ...
           median(absErr_UnLinkSelc_itr(absErr_UnLinkSelc_itr ~= -1))];
   
medianAllNormAbsErrs = [median(absErr_LinkEpi_clusterComb_itr(absErr_LinkEpi_clusterComb_itr ~= -1)) ...
           median(absErr_LinkEpi_clusterCombOnlyWell_itr(absErr_LinkEpi_clusterCombOnlyWell_itr ~= -1)) ...
           median(absErr_LinkEpi_clusterCombOnlyAmb_itr(absErr_LinkEpi_clusterCombOnlyAmb_itr ~= -1)) ...
           median(absErr_LinkEpi_cluster_itr(absErr_LinkEpi_cluster_itr ~= -1)) ...
           median(absErr_LinkEpi_clusterOnlyWell_itr(absErr_LinkEpi_clusterOnlyWell_itr ~= -1)) ...
           median(absErr_LinkEpi_clusterOnlyAmb_itr(absErr_LinkEpi_clusterOnlyAmb_itr ~= -1)) ...
           median(absErr_LinkEpi_itr(absErr_LinkEpi_itr ~= -1)) ...
           median(absErr_LinkEpi_only_si_itr(absErr_LinkEpi_only_si_itr ~= -1)) ...
           median(absErr_LinkEpi_only_sij_itr(absErr_LinkEpi_only_sij_itr ~= -1)) ...
           median(absErr_Link_itr(absErr_Link_itr ~= -1)) ...
           median(absErr_UnLink_itr(absErr_UnLink_itr ~= -1))];
%%

disp('-------------------------------')
disp(' mean AUROC   /     AUROCSelc')
allAuc = [mean(aucLinkEpiItr(aucLinkEpiItr(:,1) ~= -1, 1)) mean(aucLinkEpiItr(aucLinkEpiItr(:,2) ~= -1, 2)) mean(aucLinkEpiSelcItr(aucLinkEpiSelcItr(:,1) ~= -1,1)) mean(aucLinkEpiSelcItr(aucLinkEpiSelcItr(:,2) ~= -1,2));
mean(aucLinkNoMuEpiItr(aucLinkNoMuEpiItr(:,1) ~= -1, 1)) mean(aucLinkNoMuEpiItr(aucLinkNoMuEpiItr(:,2) ~= -1, 2)) mean(aucLinkNoMuEpiSelcItr(aucLinkNoMuEpiSelcItr(:,1) ~= -1,1)) mean(aucLinkNoMuEpiSelcItr(aucLinkNoMuEpiSelcItr(:,2) ~= -1,2));
%mean(aucLinkEpiNoEItr) mean(aucLinkEpiNoESelcItr);
mean(aucLinkItr(aucLinkItr(:,1) ~= -1, 1)) mean(aucLinkItr(aucLinkItr(:,2) ~= -1, 2)) mean(aucLinkSelcItr(aucLinkSelcItr(:,1) ~= -1,1)) mean(aucLinkSelcItr(aucLinkSelcItr(:,2) ~= -1,2));
mean(aucUnLinkItr(aucUnLinkItr(:,1) ~= -1, 1)) mean(aucUnLinkItr(aucUnLinkItr(:,2) ~= -1, 2)) mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,1) ~= -1,1)) mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,2) ~= -1,2));
mean(aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,2) ~= -1,2)) mean(aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,2) ~= -1,2));
mean(aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,2) ~= -1,2)) mean(aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,2) ~= -1,2));
mean(aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,2) ~= -1,2));
mean(aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,2) ~= -1,2));

]
disp('-------------------------------')

if(recombination == 1)
    fileNameSave = ['PlotData_WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(numItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
else
    fileNameSave = ['PlotData_WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(numItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr1_' num2str(numItr) '.mat'];
end

disp('Saving file...')
save([dirNameFigures fileNameSave])
