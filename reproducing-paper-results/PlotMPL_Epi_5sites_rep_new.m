%%

% Last updated 03-Jan-2021
% code for 5-site simulation Figure 5

% Previous updated 31-Oct 2017
% Author: M Saqib Sohail


clc
%close all
clear all

init_ColorNames_Settings;
% set(0,'DefaultAxesFontName','CMU Serif Roman')
% set(0,'DefaultTextFontName','CMU Serif Roman')
% set(0,'DefaultAxesFontSize',15)
% set(0,'DefaultTextFontSize',15)

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

exportData2Excel = 0;
saveFigs = 0;
%Lin = 10;
thisSet = 1062001%10603%5611001%46%58%1991;%1990;%42%87%58%5%87%5;
fileNameContainingDirPath = 'dirNames.txt';
postFiltering = 1; % 0 : pre filtering (out put of AnalysisMPL_filt)
                   % 1 : post filtering
similarityThresh = 0.1; % threshold that sets allowed difference between two numbers to be still declared as same value, ie., 2 = 2.01
noiseThresh = 0.0;
getSysParam;
Tstart = 1;
Tused = 100;
ng = 100;
dT = 10;
numStrainsInInitialPop = 5; % 5, 10, 30

dataFilesWith2Regs = 1;
if(dataFilesWith2Regs == 1)
     regStr1 = '1';
     regStr2 = '1';
end

thisItrStart = 1;
thisItrEnd = 3;%5;
lastItr = 999;%1000;
numRepComb = thisItrEnd - thisItrStart + 1;
numItr = lastItr - thisItrStart + 1;

numItr_rep = numItr/numRepComb;
if(numItr_rep - round(numItr_rep) ~= 0)
   disp('Error: numITr is not divisible by numRepComb...') 
   pause
end

numParam = Lin*(Lin+1)/2;
step = 0.002;
edges = [-10 -0.3:step:0.3 10] + step/2;


%Tend = Tused + Tstart - 1;
Tend = Tused + Tstart;

actualT = T/1000;
posSelc = selVal/2/Nin;
negSelc = delSelVal/2/Nin;
lineCol = {'r.-', 'b.-', 'k.-', 'g.-'};
numTps = Tused/dT;

dirNameScriptFile = pwd;    
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    display('Error: system si not unix and not PC...')
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
dirNameFigures_rep = [dirNameFigures 'Rep' chosenSlash];
if(exist(dirNameFigures, 'dir') == 0)
    mkdir(dirNameFigures)        
end
if(exist(dirNameFigures_rep, 'dir') == 0)
    mkdir(dirNameFigures_rep)        
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
aucLinkEpiItr = -1*ones(numItr,2);
aucLinkNoMuEpiItr = -1*ones(numItr,2);
aucLinkEpiNoEItr = -1*ones(numItr,2);
aucLinkItr = -1*ones(numItr,2);
aucUnLinkItr = -1*ones(numItr,2);
aucLinkEpiSelcItr = -1*ones(numItr,2);
aucLinkNoMuEpiSelcItr = -1*ones(numItr,2);
aucLinkEpiNoESelcItr = -1*ones(numItr,2);
aucLinkSelcItr = -1*ones(numItr,2);
aucUnLinkSelcItr = -1*ones(numItr,2);

numSelcSitesItr = zeros(1, numItr);
numSelcSitesEpiItr = zeros(1, numItr);
numSelcSitesEpiItr_si = zeros(1, numItr);
numSelcSitesEpiItr_sij = zeros(1, numItr);

allEstEpiItr = zeros(numItr, numParam);
allEstNoMuEpiItr = zeros(numItr, numParam);
allEstEpiNoEItr = zeros(numItr, numParam);
allEstItr = zeros(numItr, Lin);

allSelcSites = zeros(numItr, Lin);
allSelcSitesEpi = zeros(numItr, numParam);

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


%----------rep param
selcCoeffNameAllItr = repmat(' ', numItr_rep*60, 4);
estimatesAllItr = repmat(' ', numItr_rep*60, numParam);
methodsAllItr = repmat(' ', numItr_rep*60, 9);
validAllItr = repmat(' ', numItr_rep*60, 1);

aucLinkEpiItr_rep = -1*ones(numItr_rep,2);
aucLinkNoMuEpiItr_rep = -1*ones(numItr_rep,2);
aucLinkEpiNoEItr_rep = -1*ones(numItr_rep,2);
aucLinkItr_rep = -1*ones(numItr_rep,2);
aucUnLinkItr_rep = -1*ones(numItr_rep,2);
aucLinkEpiSelcItr_rep = -1*ones(numItr_rep,2);
aucLinkNoMuEpiSelcItr_rep = -1*ones(numItr_rep,2);
aucLinkEpiNoESelcItr_rep = -1*ones(numItr_rep,2);
aucLinkSelcItr_rep = -1*ones(numItr_rep,2);
aucUnLinkSelcItr_rep = -1*ones(numItr_rep,2);

numSelcSitesItr_rep = zeros(1, numItr_rep);
numSelcSitesEpiItr_rep = zeros(1, numItr_rep);
numSelcSitesEpiItr_si_rep = zeros(1, numItr_rep);
numSelcSitesEpiItr_sij_rep = zeros(1, numItr_rep);

allEstEpiItr_rep = zeros(numItr_rep, numParam);
allEstNoMuEpiItr_rep = zeros(numItr_rep, numParam);
allEstEpiNoEItr_rep = zeros(numItr_rep, numParam);
allEstItr_rep = zeros(numItr_rep, Lin);

allSelcSites_rep = zeros(numItr_rep, Lin);
allSelcSitesEpi_rep = zeros(numItr_rep, numParam);

absErr_LinkEpi_Selc_itr_rep = -1*ones(numItr_rep,1);
absErr_LinkEpi_only_si_Selc_itr_rep = -1*ones(numItr_rep,1);
absErr_LinkEpi_only_sij_Selc_itr_rep = -1*ones(numItr_rep,1);
absErr_LinkSelc_itr_rep = -1*ones(numItr_rep,1);
absErr_UnLinkSelc_itr_rep = -1*ones(numItr_rep,1);
absErr_LinkEpi_cluster_itr_rep = -1*ones(numItr_rep,1);
absErr_LinkEpi_clusterComb_itr_rep = -1*ones(numItr_rep,1);
absErr_LinkEpi_clusterOnlyWell_itr_rep = -1*ones(numItr_rep,1);
absErr_LinkEpi_clusterOnlyAmb_itr_rep = -1*ones(numItr_rep,1);
absErr_LinkEpi_clusterCombOnlyWell_itr_rep = -1*ones(numItr_rep,1);
absErr_LinkEpi_clusterCombOnlyAmb_itr_rep = -1*ones(numItr_rep,1);

numTerms_LinkEpi_cluster_itr_rep  = -1*ones(numItr_rep,1);

xCount = 1;

for itr = 1:numItr%[1:50 52:100]%1:numItr
    itr
    selcSites = [];
    if(Tstart == 1)
        fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
    else        
        fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
    end
    load([dirNameAnalysis fileName], 'Lin', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
             'sigmaEstOutLink', 'sigmaEstOutLinkEpi', 'sigmaEstOutLinkNoMuEpi',...
             'sigmaEstOutUnLink', 'sigmaEstOutLinkEpiNoE', ...
             'q', 'q11', 'priorConst', 'timeWholeCodeMPL',...'perSiteSelctionSelc', 'selcSites', ... 'sigmaEstOutLinkSelc', 'sigmaEstOutLinkNoMuSelc','sigmaEstOutUnLinkSelc', 'sigmaEstOutUnLinkNoMuSelc', ...
             'deltaLogLikeli', 'deltaLogLikeli50', 'deltaLogLikeli100',...
             'deltaLogLikeliMPL', 'deltaLogLikeliSL', 'Hmat');

    trueClassROC = 3*ones(1, Lin);
    trueClassROC(perSiteSelction == posSelc) = 4;
    trueClassROC(perSiteSelction == negSelc) = 5;
    
    perSiteSelctionAllEpiTerms = [];
    for l = 1:Lin
        perSiteSelctionAllEpiTerms = [perSiteSelctionAllEpiTerms perSiteSelctionEpi(l,l+1:end)];
    end

    perSiteSelctionAllEpiTerms = [perSiteSelction perSiteSelctionAllEpiTerms];

    trueClassROCEpi = 3*ones(1, (Lin*(Lin + 1)/2));
    trueClassROCEpi(perSiteSelction == posSelc) = 4;
    trueClassROCEpi(perSiteSelction == negSelc) = 5;
    
    [selcSitesEpi, wellCondSites, illCondSites, cluster] = findIndColsOfHmat(Hmat*dT);

    if(postFiltering == 1)
        selcSites = findIndColsOfHmat(Hmat(1:Lin,1:Lin)*dT);
        numSelcSitesItr(itr) = sum(selcSites);
        numSelcSitesEpiItr(itr) = sum(selcSitesEpi);
        numSelcSitesEpiItr_si(itr) = sum(selcSitesEpi(1:Lin));
        numSelcSitesEpiItr_sij(itr) = sum(selcSitesEpi(Lin+1:end));
        
        perSiteSelctionSelc = perSiteSelction(selcSites);
        perSiteSelctionAllEpiTermsSelc = perSiteSelctionAllEpiTerms(selcSitesEpi);
        trueClassROCSelc = trueClassROC(selcSites);
        trueClassROCEpiSelc = trueClassROCEpi(selcSitesEpi);

        sigmaEstOutLinkEpiSelc = sigmaEstOutLinkEpi(selcSitesEpi);
        sigmaEstOutLinkNoMuEpiSelc = sigmaEstOutLinkNoMuEpi(selcSitesEpi);
        sigmaEstOutLinkEpiNoESelc = sigmaEstOutLinkEpiNoE(selcSitesEpi);
        sigmaEstOutLinkSelc = sigmaEstOutLink(selcSites);
        sigmaEstOutUnLinkSelc = sigmaEstOutUnLink(selcSites);
    end
    allSelcSites(itr,:) = selcSites;
    allSelcSitesEpi(itr,:) = selcSitesEpi;
    
    timeMPLItr(itr) = timeWholeCodeMPL(itr);

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
    [aucLinkNoMuEpiItrTemp, aucLinkNoMuEpiSelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkNoMuEpi, sigmaEstOutLinkNoMuEpiSelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);
    [aucLinkEpiNoEItrTemp, aucLinkEpiNoESelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkEpiNoE, sigmaEstOutLinkEpiNoESelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);

    % only s_i term sof Epi
    sitesEpi_only_si = logical([ones(1, Lin), zeros(1, Lin*(Lin-1)/2)]);
    selcSitesEpi_only_si = selcSitesEpi & sitesEpi_only_si;
    perSiteSelction_only_si = perSiteSelctionAllEpiTerms(sitesEpi_only_si);
    perSiteSelction_only_si_Selc = perSiteSelctionAllEpiTerms(selcSitesEpi_only_si);
    sigmaEstOutLinkEpi_only_si = sigmaEstOutLinkEpi(sitesEpi_only_si);
    sigmaEstOutLinkEpi_only_si_Selc = sigmaEstOutLinkEpi(selcSitesEpi_only_si);
 
    posOnlyItrTemp_only_si = perSiteSelction_only_si > neutralLimit; %posOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
    negOnlyItrTemp_only_si = perSiteSelction_only_si < -neutralLimit; %negOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
    posOnlySelcItrTemp_only_si = perSiteSelction_only_si_Selc > neutralLimit;
    negOnlySelcItrTemp_only_si = perSiteSelction_only_si_Selc < -neutralLimit;
    [aucLinkEpiItrTemp_only_si, aucLinkEpiSelcItrTemp_only_si] = calcAUROC(perSiteSelction_only_si, perSiteSelction_only_si_Selc, sigmaEstOutLinkEpi_only_si, sigmaEstOutLinkEpi_only_si_Selc, posOnlyItrTemp_only_si, negOnlyItrTemp_only_si, posOnlySelcItrTemp_only_si, negOnlySelcItrTemp_only_si);
 
    % only s_ij term sof Epi
    sitesEpi_only_sij = logical([zeros(1, Lin), ones(1, Lin*(Lin-1)/2)]);
    selcSitesEpi_only_sij = selcSitesEpi & sitesEpi_only_sij;
    perSiteSelction_only_sij = perSiteSelctionAllEpiTerms(sitesEpi_only_sij);
    perSiteSelction_only_sij_Selc = perSiteSelctionAllEpiTerms(selcSitesEpi_only_sij);
    sigmaEstOutLinkEpi_only_sij = sigmaEstOutLinkEpi(sitesEpi_only_sij);
    sigmaEstOutLinkEpi_only_sij_Selc = sigmaEstOutLinkEpi(selcSitesEpi_only_sij);
 
    posOnlyItrTemp_only_sij = perSiteSelction_only_sij > neutralLimit; %posOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
    negOnlyItrTemp_only_sij = perSiteSelction_only_sij < -neutralLimit; %negOnlyItrTemp_AllEpiTerms & logical([ones(1, 5), zeros(1, 10)]);
    posOnlySelcItrTemp_only_sij = perSiteSelction_only_sij_Selc > neutralLimit;
    negOnlySelcItrTemp_only_sij = perSiteSelction_only_sij_Selc < -neutralLimit;
    [aucLinkEpiItrTemp_only_sij, aucLinkEpiSelcItrTemp_only_sij] = calcAUROC(perSiteSelction_only_sij, perSiteSelction_only_sij_Selc, sigmaEstOutLinkEpi_only_sij, sigmaEstOutLinkEpi_only_sij_Selc, posOnlyItrTemp_only_sij, negOnlyItrTemp_only_sij, posOnlySelcItrTemp_only_sij, negOnlySelcItrTemp_only_sij);


    
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

    aucLinkEpi_only_si_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_si;
    aucLinkEpi_only_sij_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_sij;


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
    
    %jj = jj + 1;
%-------------------------------------------------------------------------------------------    
    % 1.2 Load data (combining replicates)
    if(rem(itr, thisItrEnd) == 0)    
        
        if(Tstart == 1)
            fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff'  '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str((xCount-1)*numRepComb +thisItrStart) '_' num2str(xCount*numRepComb) '.mat'];
        else        
            fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str((xCount-1)*numRepComb +thisItrStart) '_' num2str(xCount*numRepComb) '.mat'];
        end

        load([dirNameAnalysis fileName], 'Lin', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
             'sigmaEstOutLinkRep', 'sigmaEstOutLinkEpiRep', 'sigmaEstOutLinkNoMuEpiRep',...
             'sigmaEstOutUnLinkRep', 'sigmaEstOutLinkEpiNoERep', 'HmatDTAll');

         
        sigmaEstOutLink_rep = sigmaEstOutLinkRep;
        sigmaEstOutLinkEpi_rep = sigmaEstOutLinkEpiRep;
        sigmaEstOutLinkNoMuEpi_rep = sigmaEstOutLinkNoMuEpiRep;
        sigmaEstOutUnLink_rep = sigmaEstOutUnLinkRep;
        sigmaEstOutLinkEpiNoE_rep = sigmaEstOutLinkEpiNoERep;
         
         
        %==================================================================
        
        [selcSitesEpi_rep, wellCondSites_rep, illCondSites_rep, cluster_rep] = findIndColsOfHmat(HmatDTAll);
        %selcSites_rep = sum(allSelcSites((xCount-1)*numRepComb + 1:xCount*numRepComb,:)) > 0;
        selcSites_rep = findIndColsOfHmat(HmatDTAll(1:Lin,1:Lin));
        
        numSelcSitesItr_rep(xCount) = sum(selcSites_rep);
        numSelcSitesEpiItr_rep(xCount) = sum(selcSitesEpi_rep);
        numSelcSitesEpiItr_si_rep(xCount) = sum(selcSitesEpi_rep(1:Lin));
        numSelcSitesEpiItr_sij_rep(xCount) = sum(selcSitesEpi_rep(Lin+1:end));




        perSiteSelctionSelc_rep = perSiteSelction(selcSites_rep);
        perSiteSelctionAllEpiTermsSelc_rep = perSiteSelctionAllEpiTerms(selcSitesEpi_rep);

        sigmaEstOutLinkEpiSelc_rep = sigmaEstOutLinkEpi(selcSitesEpi_rep);
        sigmaEstOutLinkNoMuEpiSelc_rep = sigmaEstOutLinkNoMuEpi(selcSitesEpi_rep);
        sigmaEstOutLinkEpiNoESelc_rep = sigmaEstOutLinkEpiNoE(selcSitesEpi_rep);
        sigmaEstOutLinkSelc_rep = sigmaEstOutLink(selcSites_rep);
        sigmaEstOutUnLinkSelc_rep = sigmaEstOutUnLink(selcSites_rep);

        
        
        
        [aucLinkItrTemp_rep, aucLinkSelcItrTemp_rep, aucUnLinkItrTemp_rep, aucUnLinkSelcItrTemp_rep, ...
        aucLinkEpiItrTemp_rep, aucLinkEpiSelcItrTemp_rep, aucLinkNoMuEpiItrTemp_rep, aucLinkNoMuEpiSelcItrTemp_rep, ...    
        aucLinkEpiNoEItrTemp_rep, aucLinkEpiNoESelcItrTemp_rep, aucLinkEpiItrTemp_only_si_rep, aucLinkEpiSelcItrTemp_only_si_rep, ...
        aucLinkEpiItrTemp_only_sij_rep, aucLinkEpiSelcItrTemp_only_sij_rep] = calcEpiRef_AUROC_new(perSiteSelction, perSiteSelctionSelc_rep, ...
                            perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc_rep, selcSitesEpi_rep, ...
                            sigmaEstOutLink_rep, sigmaEstOutLinkSelc_rep, sigmaEstOutUnLink_rep, sigmaEstOutUnLinkSelc_rep, ...
                            sigmaEstOutLinkEpi_rep, sigmaEstOutLinkEpiSelc_rep, sigmaEstOutLinkNoMuEpi_rep, sigmaEstOutLinkNoMuEpiSelc_rep, ...
                            sigmaEstOutLinkEpiNoE_rep, sigmaEstOutLinkEpiNoESelc_rep, neutralLimit);
                    
                        
                        
        if(sum(selcSites_rep) > 0)
            if(sum(abs(perSiteSelction(selcSites_rep))) < zeroThresh)
                absErr_LinkSelc_itr_rep(xCount) = sum(abs(perSiteSelction(selcSites_rep) - sigmaEstOutLinkSelc_rep'))/(sum(abs(perSiteSelction(selcSites_rep))) + naeReg);%/length(perSiteSelction(selcSites_rep));
                absErr_UnLinkSelc_itr_rep(xCount) = sum(abs(perSiteSelction(selcSites_rep) - sigmaEstOutUnLinkSelc_rep'))/(sum(abs(perSiteSelction(selcSites_rep))) + naeReg);%/length(perSiteSelction(selcSites_rep));
            else
                absErr_LinkSelc_itr_rep(xCount) = sum(abs(perSiteSelction(selcSites_rep) - sigmaEstOutLinkSelc_rep'))/sum(abs(perSiteSelction(selcSites_rep)));%/length(perSiteSelction(selcSites_rep));
                absErr_UnLinkSelc_itr_rep(xCount) = sum(abs(perSiteSelction(selcSites_rep) - sigmaEstOutUnLinkSelc_rep'))/sum(abs(perSiteSelction(selcSites_rep)));%/length(perSiteSelction(selcSites_rep));
            end
        end
        
       
            % Calculate ABSOLUTE ERROR metrics
    
    
    
    
        if(sum(abs(perSiteSelction)) < zeroThresh)
            absErr_Link_itr_rep(xCount) = sum(abs(perSiteSelction - sigmaEstOutLink_rep'))/(sum(abs(perSiteSelction)) + naeReg);%/length(perSiteSelction);
            absErr_UnLink_itr_rep(xCount) = sum(abs(perSiteSelction - sigmaEstOutUnLink_rep'))/(sum(abs(perSiteSelction)) + naeReg);%/length(perSiteSelction);
        else
            absErr_Link_itr_rep(xCount) = sum(abs(perSiteSelction - sigmaEstOutLink_rep'))/sum(abs(perSiteSelction));%/length(perSiteSelction);
            absErr_UnLink_itr_rep(xCount) = sum(abs(perSiteSelction - sigmaEstOutUnLink_rep'))/sum(abs(perSiteSelction));%/length(perSiteSelction);
        end

        perSiteSelctionAllEpiTerms_cluster_rep = perSiteSelctionAllEpiTerms(wellCondSites_rep);
        sigmaEstOutLinkEpi_cluster_rep = sigmaEstOutLinkEpi_rep(wellCondSites_rep)';
        perSiteSelctionAllEpiTerms_clusterComb_rep = perSiteSelctionAllEpiTerms_cluster_rep;
        sigmaEstOutLinkEpi_clusterComb_rep = sigmaEstOutLinkEpi_cluster_rep;
        indTemp201_rep = [];
        num_clusterWellOnly_rep = length(sigmaEstOutLinkEpi_cluster_rep);
        for wt = 1:length(cluster_rep)
            indTemp101_rep = cluster_rep{wt};
            indTemp201_rep = [indTemp201_rep indTemp101_rep];
            perSiteSelctionAllEpiTerms_cluster_rep = [perSiteSelctionAllEpiTerms_cluster_rep sum(perSiteSelctionAllEpiTerms(indTemp101_rep))];
            sigmaEstOutLinkEpi_cluster_rep = [sigmaEstOutLinkEpi_cluster_rep sum(sigmaEstOutLinkEpi_rep(indTemp101_rep))];
        end
        clustInd_rep = zeros(1, length(sigmaEstOutLinkEpi_cluster_rep));
        temp65_rep = 1:num_clusterWellOnly_rep;
        clustInd_rep(temp65_rep) = 1;
        clustInd_rep = logical(clustInd_rep);



        perSiteSelctionAllEpiTerms_cluster_rep = round(100000*perSiteSelctionAllEpiTerms_cluster_rep)/100000;
        if(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep)) < zeroThresh)
            absErr_LinkEpi_cluster_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep - sigmaEstOutLinkEpi_cluster_rep))/sum(abs(perSiteSelctionAllEpiTerms_cluster_rep) + naeReg);%/length(perSiteSelctionAllEpiTerms_cluster_rep);
        else
            absErr_LinkEpi_cluster_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep - sigmaEstOutLinkEpi_cluster_rep))/sum(abs(perSiteSelctionAllEpiTerms_cluster_rep));%/length(perSiteSelctionAllEpiTerms_cluster_rep);
        end
        numTerms_LinkEpi_cluster_itr_rep(xCount) = length(perSiteSelctionAllEpiTerms_cluster_rep);

        cluster_detail_cell_rep{itr,1} = perSiteSelctionAllEpiTerms_cluster_rep;
        cluster_detail_cell_rep{itr,2} = sigmaEstOutLinkEpi_cluster_rep;

        if(sum(clustInd_rep) > 0)
            if(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep))) < zeroThresh)
                absErr_LinkEpi_clusterOnlyWell_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep) - sigmaEstOutLinkEpi_cluster_rep(clustInd_rep)))/(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep))) + naeReg);%/sum(clustInd_rep);
            else
                absErr_LinkEpi_clusterOnlyWell_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep) - sigmaEstOutLinkEpi_cluster_rep(clustInd_rep)))/sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep)));%/sum(clustInd_rep);
            end
        end

        if(sum(~clustInd_rep) > 0)
            if(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep))) < zeroThresh)
                absErr_LinkEpi_clusterOnlyAmb_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep) - sigmaEstOutLinkEpi_cluster_rep(~clustInd_rep)))/(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep))) + naeReg);%/sum(~clustInd_rep);
            else
                absErr_LinkEpi_clusterOnlyAmb_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep) - sigmaEstOutLinkEpi_cluster_rep(~clustInd_rep)))/sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep)));%/sum(~clustInd_rep);
            end
        end

        if(isinf(absErr_LinkEpi_cluster_itr_rep(xCount)))
            pause
        end

        if(absErr_LinkEpi_cluster_itr_rep(xCount) > 5)
            %pause
        end

        % cluster     : well cond sites unique and ill cond sites divided into 
        %               amb clusters, each cluster as 1 selCoeff
        % clusterComb : all amb cluster bunched in 1 selCoeff

        
        allSelcSites_rep(xCount,:) = selcSites_rep;
        allSelcSitesEpi_rep(xCount,:) = selcSitesEpi_rep;
         
        aucLinkEpiItr_rep(xCount,:) = aucLinkEpiItrTemp_rep;
        aucLinkNoMuEpiItr_rep(xCount,:) = aucLinkNoMuEpiItrTemp_rep;
        aucLinkEpiNoEItr_rep(xCount,:) = aucLinkEpiNoEItrTemp_rep;
        aucLinkItr_rep(xCount,:) = aucLinkItrTemp_rep;
        aucUnLinkItr_rep(xCount,:) = aucUnLinkItrTemp_rep;

        % calculations for filtered
        aucLinkEpiSelcItr_rep(xCount,:) = aucLinkEpiSelcItrTemp_rep;
        aucLinkNoMuEpiSelcItr_rep(xCount,:) = aucLinkNoMuEpiSelcItrTemp_rep;
        aucLinkEpiNoESelcItr_rep(xCount,:) = aucLinkEpiNoESelcItrTemp_rep;
        aucLinkSelcItr_rep(xCount,:) = aucLinkSelcItrTemp_rep;
        aucUnLinkSelcItr_rep(xCount,:) = aucUnLinkSelcItrTemp_rep;

        aucLinkEpi_only_si_Itr_rep(xCount,:) = aucLinkEpiItrTemp_only_si_rep;
        aucLinkEpi_only_sij_Itr_rep(xCount,:) = aucLinkEpiItrTemp_only_sij_rep;

        aucLinkEpi_only_si_SelcItr_rep(xCount,:) = aucLinkEpiSelcItrTemp_only_si_rep;
        aucLinkEpi_only_sij_SelcItr_rep(xCount,:) = aucLinkEpiSelcItrTemp_only_sij_rep;

        condHmatWellItr_rep(itr) = cond(HmatDTAll(wellCondSites_rep,wellCondSites_rep));
        
        allEstEpiItr_rep(xCount,:) = sigmaEstOutLinkEpi_rep';
        allEstNoMuEpiItr_rep(xCount,:) = sigmaEstOutLinkNoMuEpi_rep';
        allEstEpiNoEItr_rep(xCount,:) = sigmaEstOutLinkEpiNoE_rep';
        allEstItr_rep(xCount,:) = sigmaEstOutLink_rep';
        
        % prep data from export to R
        selcCoeffNameAllItr_rep((xCount-1)*60 + 1:xCount*60,:) = repmat(selcCoeffNameColumnOneItr, 4, 1);

        temp15 = num2str([sigmaEstOutLinkEpi_rep; sigmaEstOutLinkNoMuEpi_rep]);

        temp20 = num2str(sigmaEstOutLink_rep); 
        temp21 = repmat(' ', 10, size(temp20, 2));
        temp25 = [temp20; temp21];

        temp30 = num2str(sigmaEstOutUnLink_rep); 
        temp31 = repmat(' ', 10, size(temp30, 2));
        temp35 = [temp30; temp31];

%         %[(xCount-1)*60+1  xCount*60-30  (xCount-1)*60+1+30  xCount*60-15 (xCount-1)*60+1+45   xCount*60
%         estimatesAllItr_rep((xCount-1)*60+1:xCount*60-30,1:size(temp15,2)) = temp15;
%         estimatesAllItr_rep((xCount-1)*60+1+30:xCount*60-15,1:size(temp25,2)) = temp25;
%         estimatesAllItr_rep((xCount-1)*60+1+45:xCount*60,1:size(temp35,2)) = temp35;
% 
%         methodsAllItr_rep((xCount-1)*60 + 1:xCount*60,:) = methodCol;
% 
%         validAllItr_rep((xCount-1)*60+1:xCount*60-30,1) = num2str([selcSitesEpi_rep selcSitesEpi_rep]');
%         validAllItr_rep((xCount-1)*60+1+30:xCount*60-15,1) = [num2str(selcSites_rep'); repmat(' ', 10, 1)];
%         validAllItr_rep((xCount-1)*60+1+45:xCount*60,1) = [num2str(selcSites_rep'); repmat(' ', 10, 1)];

        xCount = xCount + 1;
    end
    
end

%% estimation performance

% % only positive sites
posOnly = perSiteSelction_AllItr == max(perSiteSelction_AllItr);
negOnly = perSiteSelction_AllItr == min(perSiteSelction_AllItr);
neutOnly = perSiteSelction_AllItr == 0;
%

%%
disp('-------------------------------')
disp(' mean AUROC   /     AUROCSelc')
allAuc = [mean(aucLinkEpiItr(aucLinkEpiItr(:,1) ~= -1, 1)) mean(aucLinkEpiItr(aucLinkEpiItr(:,2) ~= -1, 2)) mean(aucLinkEpiSelcItr(aucLinkEpiSelcItr(:,1) ~= -1,1)) mean(aucLinkEpiSelcItr(aucLinkEpiSelcItr(:,2) ~= -1,2));
mean(aucLinkNoMuEpiItr(aucLinkNoMuEpiItr(:,1) ~= -1, 1)) mean(aucLinkNoMuEpiItr(aucLinkNoMuEpiItr(:,2) ~= -1, 2)) mean(aucLinkNoMuEpiSelcItr(aucLinkNoMuEpiSelcItr(:,1) ~= -1,1)) mean(aucLinkNoMuEpiSelcItr(aucLinkNoMuEpiSelcItr(:,2) ~= -1,2));
%mean(aucLinkEpiNoEItr) mean(aucLinkEpiNoESelcItr);
mean(aucLinkItr(aucLinkItr(:,1) ~= -1, 1)) mean(aucLinkItr(aucLinkItr(:,2) ~= -1, 2)) mean(aucLinkSelcItr(aucLinkSelcItr(:,1) ~= -1,1)) mean(aucLinkSelcItr(aucLinkSelcItr(:,2) ~= -1,2));
mean(aucUnLinkItr(aucUnLinkItr(:,1) ~= -1, 1)) mean(aucUnLinkItr(aucUnLinkItr(:,2) ~= -1, 2)) mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,1) ~= -1,1)) mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,2) ~= -1,2));
mean(aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,2) ~= -1,2)) mean(aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,2) ~= -1,2));
mean(aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,2) ~= -1,2)) mean(aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,2) ~= -1,2))]
disp('-------------------------------')

% replicate comb
% disp(' mean AUROC   /     AUROCSelc')
% allAuc_rep = [mean(aucLinkEpiItr_rep) mean(aucLinkEpiSelcItr_rep);
% mean(aucLinkNoMuEpiItr_rep) mean(aucLinkNoMuEpiSelcItr_rep);
% %mean(aucLinkEpiNoEItr_rep) mean(aucLinkEpiNoESelcItr_rep);
% mean(aucLinkItr_rep) mean(aucLinkSelcItr_rep);
% mean(aucUnLinkItr_rep) mean(aucUnLinkSelcItr_rep)]


% replicate comb
disp('-------------------------------')
disp(' mean AUROC   /     AUROCSelc')
allAuc_rep = [mean(aucLinkEpiItr_rep(aucLinkEpiItr_rep(:,1) ~= -1, 1)) mean(aucLinkEpiItr_rep(aucLinkEpiItr_rep(:,2) ~= -1, 2)) mean(aucLinkEpiSelcItr_rep(aucLinkEpiSelcItr_rep(:,1) ~= -1,1)) mean(aucLinkEpiSelcItr_rep(aucLinkEpiSelcItr_rep(:,2) ~= -1,2));
mean(aucLinkNoMuEpiItr_rep(aucLinkNoMuEpiItr_rep(:,1) ~= -1, 1)) mean(aucLinkNoMuEpiItr_rep(aucLinkNoMuEpiItr_rep(:,2) ~= -1, 2)) mean(aucLinkNoMuEpiSelcItr_rep(aucLinkNoMuEpiSelcItr_rep(:,1) ~= -1,1)) mean(aucLinkNoMuEpiSelcItr_rep(aucLinkNoMuEpiSelcItr_rep(:,2) ~= -1,2));
%mean(aucLinkEpiNoEItr_rep) mean(aucLinkEpiNoESelcItr_rep);
mean(aucLinkItr_rep(aucLinkItr_rep(:,1) ~= -1, 1)) mean(aucLinkItr_rep(aucLinkItr_rep(:,2) ~= -1, 2)) mean(aucLinkSelcItr_rep(aucLinkSelcItr_rep(:,1) ~= -1,1)) mean(aucLinkSelcItr_rep(aucLinkSelcItr_rep(:,2) ~= -1,2));
mean(aucUnLinkItr_rep(aucUnLinkItr_rep(:,1) ~= -1, 1)) mean(aucUnLinkItr_rep(aucUnLinkItr_rep(:,2) ~= -1, 2)) mean(aucUnLinkSelcItr_rep(aucUnLinkSelcItr_rep(:,1) ~= -1,1)) mean(aucUnLinkSelcItr_rep(aucUnLinkSelcItr_rep(:,2) ~= -1,2));
mean(aucLinkEpi_only_si_Itr_rep(aucLinkEpi_only_si_Itr_rep(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_Itr_rep(aucLinkEpi_only_si_Itr_rep(:,2) ~= -1,2)) mean(aucLinkEpi_only_si_SelcItr_rep(aucLinkEpi_only_si_SelcItr_rep(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_SelcItr_rep(aucLinkEpi_only_si_SelcItr_rep(:,2) ~= -1,2));
mean(aucLinkEpi_only_sij_Itr_rep(aucLinkEpi_only_sij_Itr_rep(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_Itr_rep(aucLinkEpi_only_sij_Itr_rep(:,2) ~= -1,2)) mean(aucLinkEpi_only_sij_SelcItr_rep(aucLinkEpi_only_sij_SelcItr_rep(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_SelcItr_rep(aucLinkEpi_only_sij_SelcItr_rep(:,2) ~= -1,2))]
disp('-------------------------------')




%%
%%
% bar plot of MPLE and MPL
%----------------------------
for k = 1:2
    fig1 = figure
    sub1 = subplot(1,2,1)
    pause(1)
    
    if(k == 1)
        %barDataToPlot = allAuc([3 4 5 6 1],[1 2])';% unfiltered   %[aucLinkItr ;aucLinkNoMuItr; aucUnLinkItr; aucUnLinkNoMuItr]';
        barDataToPlot = allAuc([3 4 5 6],[1 2])';% unfiltered   %[aucLinkItr ;aucLinkNoMuItr; aucUnLinkItr; aucUnLinkNoMuItr]';
        titleStr = ' Auc';
        figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_AUROC_allSites' '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
    else
        %barDataToPlot = allAuc([3 4 5 6 1],[3 4])';% Selc
        barDataToPlot = allAuc_rep([3 4 5 6],[1 2])';% Selc
        titleStr = 'Rep';
        figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_AUROC_selcSites' '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
    end

    bar(barDataToPlot)
    %bar(barDataToPlot(:,[1 4]))
    %[map,num,typ] = brewermap(4, 'Paired');
    map6Paired = brewermap(8,'Paired');
    mapOranges = brewermap(8,'Oranges');
    mapBlues = brewermap(8,'Blues');
    myColorMap = color_scheme_npg([2:8],:);
    %color_scheme_npg
    %colormap([mapOranges(5,:); mapOranges(3,:); mapBlues(5,:); mapBlues(3,:)])
    colormap([mapOranges(5,:); mapBlues(5,:); mapBlues(3,:)])
    colormap(myColorMap);
    % map = [0.651 0.8078 0.8902;
    %           0.1216 0.4706 0.7059;
    %           0.6980 0.8745 0.5412;
    %           0.2000 0.6275 0.1725];
    % 
    % colormap( map)

    LabelBen = { 'Estimator (beneficial)'};
    LabelDele = {'Estimator (deleterious)'};
    set(sub1, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...
        'XTickLabel', LabelBen, ...'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...
        'YTick'       , 0:0.1:1, ...'XTick'       , 0.5:0.5:1.5, ...
      'LineWidth', 1)
    ylabel('AUROC')
    % xlabel('            ')
    axis([0.5 1.5 0.5 1])
    title(titleStr)
    %leg = legend('Multi-locus diffusion', 'Multi-locus diffusion without \mu component', 'Multi-locus diffusion without linkage component', 'Multi-locus diffusion without \mu and linkage components', 'location', 'NorthOutSide');
    %leg = legend('Multi-locus diffusion with \mu component', 'Multi-locus diffusion without \mu component', 'Single-locus diffusion with \mu component', 'Single-locus diffusion without \mu component', 'location', 'NorthOutSide');
    %leg = legend('MPL with \mu component', 'MPL without \mu component', 'Single-locus diffusion with \mu component', 'Single-locus diffusion without \mu component', 'location', 'NorthOutSide');
    leg = legend('MPL', 'SL', 'MPLE-s_i', 'MPLE-s_{ij}', 'MPLE-All', 'location', 'NorthOutSide');
    %leg = legend('MPLE', 'MPL', 'SL', 'location', 'NorthEast');

    set(leg,'color','none');
    set(leg, 'Edgecolor','none');
    % title(['Selection: ' num2str(selVal/2/Nin)])
    pause(1)
    %pos1 = get(gca, 'Position')

    set(leg, 'Visible', 'off')
    subplot(1,2,2)
    %subplot(1,4,4)

    bar(barDataToPlot)
    %bar(barDataToPlot(:,[1 4]));

    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...
      'XTickLabel', LabelDele, ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...
      'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
      'LineWidth', 1)
    % ylabel('AUROC')
    axis([1.5 2.5 0.5 1])
    title(titleStr)
    %leg = legend('Multi-locus diffusion', 'Multi-locus diffusion without \mu component', 'Multi-locus diffusion without linkage component', 'Multi-locus diffusion without \mu and linkage components', 'location', 'NorthOutSide');
    %leg = legend('Multi-locus diffusion with \mu component', 'Multi-locus diffusion without \mu component', 'Single-locus diffusion with \mu component', 'Single-locus diffusion without \mu component', 'location', 'NorthOutSide');
    %leg = legend('MPL with \mu component', 'MPL without \mu component', 'Single-locus diffusion with \mu component', 'Single-locus diffusion without \mu component', 'location', 'NorthOutSide');
    %leg = legend('MPLE', 'MPLE without \mu', 'MPL', 'SL', 'location', 'NorthOutSide');
    leg = legend('MPL', 'SL', 'MPLE-s_i', 'MPLE-s_{ij}', 'MPLE-All', 'location', 'NorthOutSide');
    set(leg,'color','none');
    set(leg, 'Edgecolor','none');
    
    pause(1)
    % title(' ')

end





%%
% bar plot of ONLY MPLE 
%----------------------------
fig5 = figure
pause(1)

%barDataToPlot = allAuc(:,[1 2])';% unfiltered   %[aucLinkItr ;aucLinkNoMuItr; aucUnLinkItr; aucUnLinkNoMuItr]';
barDataToPlot = allAuc(1,[3 4])';% Selc
barDataToPlot = [barDataToPlot barDataToPlot];
bar(barDataToPlot)

%myColorMap = color_scheme_npg([4 2 8],:);
%color_scheme_npg
%colormap([mapOranges(5,:); mapOranges(3,:); mapBlues(5,:); mapBlues(3,:)])
%colormap([mapOranges(5,:); mapBlues(5,:); mapBlues(3,:)])
alpha1 = 0.6;
myColorMap2(1,:) = color_scheme_npg(4,:)*(alpha1) + [1 1 1]*(1 - alpha1);
myColorMap2(2,:) = color_scheme_npg(2,:)*(alpha1) + [1 1 1]*(1 - alpha1);
myColorMap2(3,:) = color_scheme_npg(8,:)*(alpha1) + [1 1 1]*(1 - alpha1);
colormap(myColorMap2);




% map = [0.651 0.8078 0.8902;
%           0.1216 0.4706 0.7059;
%           0.6980 0.8745 0.5412;
%           0.2000 0.6275 0.1725];
% 
% colormap( map)


LabelBen = {'Performance'};
set(gca, ...
    'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...'XTick', 1:sum(indSelected_logical), ...
  'XTickLabel', LabelBen, ...'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...
    'YTick'       , 0:0.1:1)%, ...'XTick'       , 0.5:0.5:1.5, ...
  %'LineWidth', 1)
ylabel('AUROC')
% xlabel('            ')
axis([0.5 1.5 0.5 1])
leg = legend('Beneficial', 'Deleterious', 'location', 'NorthOutSide');

set(leg,'color','none');
set(leg, 'Edgecolor','none');

% if(saveFigs == 1)
%     %figname = ['/local/staff/ee/mssohail/Matlab Codes/H3N2 evolution/Fiftysites_3classPotts_Set' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_beneDele_all4methods_barChart'];
%     figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_beneDele_OnlyMPLE_barChart'];
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 4.5 7.2])%[0 0 19 6])% ,[0 0 8 6])
%     set(gcf, 'renderer', 'painters');
%     print(figname, '-dpng','-r400')
%     %print(figname, '-depsc')
% end

%%
% make chord plots

% myColorMap = color_scheme_npg([4 2 8],:);
% 
% alpha1 = 0.6;
% 
% myColorMap(1,:) = myColorMap(1,:)*(alpha1) + [1 1 1]*(1 - alpha1);
% myColorMap(2,:) = myColorMap(2,:)*(alpha1) + [1 1 1]*(1 - alpha1);
% myColorMap(3,:) = myColorMap(3,:)*(alpha1) + [1 1 1]*(1 - alpha1);

myColorMap5 = [myColorMap2; 0 0 0];
% Create custom node labels
myLabel = cell(Lin);
for i = 1:5
  myLabel{i} = ['s_' num2str(i)];%['L_' num2str(i)];
end

perSiteSelctionEpiRd = round(perSiteSelctionEpi*1000)/1000;
figure
circularGraph_noButtons(perSiteSelctionEpiRd*100,'Colormap',myColorMap5,'Label',myLabel);

% if(saveFigs == 1)
%     %figname = ['/local/staff/ee/mssohail/Matlab Codes/H3N2 evolution/Fiftysites_3classPotts_Set' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_beneDele_all4methods_barChart'];
%     figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_chordPlot_GT'];
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 4 4])% ,[0 0 8 6])
%     %set(gcf, 'renderer', 'painters');
%     print(figname, '-dpng','-r400')
%     %print(figname, '-depsc')
% end

%%
%selCoeffBoxPlot(allEstEpiItr, allEstItr, perSiteSelctionEpi, Lin)
allEstEpiForBoxPlot_rep = [];
allClassEpiIndicatorForBoxPlot_rep = [];
allEstForBoxPlot_rep = [];
allClassIndicatorForBoxPlot_rep = [];
for l2 = 1:15
    thisSelcSiteEpiLogicInd_rep = logical(allSelcSitesEpi_rep(:,l2));
    allEstEpiForBoxPlot_rep = [allEstEpiForBoxPlot_rep ;allEstEpiItr_rep(thisSelcSiteEpiLogicInd_rep,l2)];
    allClassEpiIndicatorForBoxPlot_rep = [allClassEpiIndicatorForBoxPlot_rep ; repmat(l2,  sum(thisSelcSiteEpiLogicInd_rep), 1)]; 
    if(l2 < 6)
        thisSelcSiteLogicInd_rep = logical(allSelcSites_rep(:,l2));
        allEstForBoxPlot_rep = [allEstForBoxPlot_rep ;allEstItr_rep(thisSelcSiteLogicInd_rep,l2)];
        allClassIndicatorForBoxPlot_rep = [allClassIndicatorForBoxPlot_rep ; repmat(l2,  sum(thisSelcSiteLogicInd_rep), 1)]; 
    end
end

color_scheme = [repmat(color_scheme_npg(3,:), 5, 1); repmat(color_scheme_npg(5,:), 10, 1)];
only1plot = 0;
selCoeffBoxPlot_new(allEstEpiForBoxPlot_rep, allClassEpiIndicatorForBoxPlot_rep, allEstForBoxPlot_rep, allClassIndicatorForBoxPlot_rep, perSiteSelctionEpi, Lin, color_scheme, only1plot)

if(saveFigs == 1)
    figname = ['FiveSites_Rep'  num2str(numRepComb) '_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_boxPlot' '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 12])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print([dirNameFigures_rep figname], '-dpng','-r400')
end

only1plot = 1;
selCoeffBoxPlot_new(allEstEpiForBoxPlot_rep, allClassEpiIndicatorForBoxPlot_rep, allEstForBoxPlot_rep, allClassIndicatorForBoxPlot_rep, perSiteSelctionEpi, Lin, color_scheme, only1plot)

if(saveFigs == 1)
    figname = ['FiveSites_Rep'  num2str(numRepComb) '_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_boxPlot_onlyMPLE' '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 4])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print([dirNameFigures_rep figname], '-dpng','-r400')
end

%%
estSiteSelctionEpi = zeros(Lin, Lin);
estSiteSelctionEpi_rep = zeros(Lin, Lin);
for l2 = 1:15
    thisSelcSiteLogicInd = logical(ones(numItr,1));%logical(allSelcSitesEpi(:,l2));
    thisSelcSiteLogicInd_rep = logical(ones(numItr_rep,1));%logical(allSelcSitesEpi_rep(:,l2));
    if(l2 <= 5)
        estSiteSelctionEpi(l2,l2) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        estSiteSelctionEpi_rep(l2,l2) = mean(allEstEpiItr_rep(thisSelcSiteLogicInd_rep,l2));
    elseif(l2 > 5 && l2 <= 9)
        estSiteSelctionEpi(1,l2-4) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        estSiteSelctionEpi_rep(1,l2-4) = mean(allEstEpiItr_rep(thisSelcSiteLogicInd_rep,l2));
    elseif(l2 > 9 && l2 <= 12)
        estSiteSelctionEpi(2,l2-7) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        estSiteSelctionEpi_rep(2,l2-7) = mean(allEstEpiItr_rep(thisSelcSiteLogicInd_rep,l2));
    elseif(l2 > 12 && l2 <= 14)
        estSiteSelctionEpi(3,l2-9) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        estSiteSelctionEpi_rep(3,l2-9) = mean(allEstEpiItr_rep(thisSelcSiteLogicInd_rep,l2));
    elseif(l2 > 14)
        estSiteSelctionEpi(4,l2-10) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        estSiteSelctionEpi_rep(4,l2-10) = mean(allEstEpiItr_rep(thisSelcSiteLogicInd_rep,l2));
    end
end
estSiteSelctionEpi(isnan(estSiteSelctionEpi)) = 0;
estSiteSelctionEpiRd = round(estSiteSelctionEpi*10000)/10000;

estSiteSelctionEpi_rep(isnan(estSiteSelctionEpi_rep)) = 0;
estSiteSelctionEpiRd_rep = round(estSiteSelctionEpi_rep*10000)/10000;

figure
circularGraph_noButtons(estSiteSelctionEpiRd*100,'Colormap',myColorMap5,'Label',myLabel);

% if(saveFigs == 1)
%     figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_chordPlot_estEpiAve'];
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 6 6])%,[0 0 8 6.45])% ,[0 0 8 6])
%     %set(gcf, 'renderer', 'painters');
%     print([dirNameFigures_rep figname], '-dpng','-r400')
% end

figure
circularGraph_noButtons(estSiteSelctionEpiRd_rep*100,'Colormap',myColorMap5,'Label',myLabel);

if(saveFigs == 1)
    figname = ['FiveSites_Rep' num2str(numRepComb) '_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_chordPlot_estEpiAve' '_reg1' regStr1  '_reg2' regStr2 '.png'];
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 4 4])%6 6])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print([dirNameFigures_rep figname], '-dpng','-r400')
end



%% make the chord diagram for a single itr (works only for Lin == 5)

chosenItrAll = 1:3%[1 2 3 4]%13;

for i = 1:length(chosenItrAll)
    chosenItr = chosenItrAll(i);
    estSiteSelctionEpiChosenItr = zeros(Lin, Lin);
    for l2 = 1:15
        thisSelcSiteLogicInd = allSelcSitesEpi(chosenItr,l2);
        if(thisSelcSiteLogicInd == 1)
            thisEntry = allEstEpiItr(chosenItr, l2);
        else
            thisEntry = 0;
        end
        if(l2 <= 5)
            estSiteSelctionEpiChosenItr(l2,l2) = thisEntry;
        elseif(l2 > 5 && l2 <= 9)
            estSiteSelctionEpiChosenItr(1,l2-4) = thisEntry;
        elseif(l2 > 9 && l2 <= 12)
            estSiteSelctionEpiChosenItr(2,l2-7) = thisEntry;
        elseif(l2 > 12 && l2 <= 14)
            estSiteSelctionEpiChosenItr(3,l2-9) = thisEntry;
        elseif(l2 > 14)
            estSiteSelctionEpiChosenItr(4,l2-10) = thisEntry;
        end
    end

    estSiteSelctionEpiChosenItrRd = round(estSiteSelctionEpiChosenItr*1000)/1000;
    figure
    circularGraph_noButtons(estSiteSelctionEpiChosenItrRd*100,'Colormap',myColorMap5,'Label',myLabel);

%     if(saveFigs == 1)
%         figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_chordPlot_estEpi_itr' num2str(chosenItr)];
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 3 3])%,[0 0 8 6.45])% ,[0 0 8 6])
%         %set(gcf, 'renderer', 'painters');
%         print([dirNameFigure figname], '-dpng','-r400')
%     end
end


%% (works only for Lin == 5)
chosenItrAll_rep = 1:2%13;

for i = 1:length(chosenItrAll_rep)
    chosenItr_rep = chosenItrAll_rep(i);
    estSiteSelctionEpiChosenItr_rep = zeros(Lin, Lin);
    for l2 = 1:15
        thisSelcSiteLogicInd_rep = allSelcSitesEpi_rep(chosenItr_rep,l2);
        if(thisSelcSiteLogicInd_rep == 1)
            thisEntry = allEstEpiItr_rep(chosenItr_rep, l2);
        else
            thisEntry = 0;
        end
        if(l2 <= 5)
            estSiteSelctionEpiChosenItr_rep(l2,l2) = thisEntry;
        elseif(l2 > 5 && l2 <= 9)
            estSiteSelctionEpiChosenItr_rep(1,l2-4) = thisEntry;
        elseif(l2 > 9 && l2 <= 12)
            estSiteSelctionEpiChosenItr_rep(2,l2-7) = thisEntry;
        elseif(l2 > 12 && l2 <= 14)
            estSiteSelctionEpiChosenItr_rep(3,l2-9) = thisEntry;
        elseif(l2 > 14)
            estSiteSelctionEpiChosenItr_rep(4,l2-10) = thisEntry;
        end
    end

    estSiteSelctionEpiChosenItrRd_rep = round(estSiteSelctionEpiChosenItr_rep*1000)/1000;
    figure
    circularGraph_noButtons(estSiteSelctionEpiChosenItrRd_rep*100,'Colormap',myColorMap5,'Label',myLabel);

    if(saveFigs == 1)
        figname = ['FiveSites_Rep' num2str(numRepComb) '_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_chordPlot_estEpi_itr' num2str(chosenItr_rep) '_reg1' regStr1  '_reg2' regStr2 '.png'];
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 3 3])%,[0 0 8 6.45])% ,[0 0 8 6])
        %set(gcf, 'renderer', 'painters');
        print([dirNameFigures_rep figname], '-dpng','-r400')
    end
end



if(exportData2Excel)
    SelCoeff = selcCoeffNameAllItr;
    Estimates = estimatesAllItr;
    Method = methodsAllItr;
    Valid = validAllItr;
    TableOfSelecEst = table(SelCoeff, Estimates, Method, Valid);
    fileName_Table = [dirNameScriptFile chosenSlash 'TableEstimates_5Sites_Set' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_numStrains' num2str(numStrainsInInitialPop) '_reg1' regStr1  '_reg2' regStr2 '.xlsx'];
    writetable(TableOfSelecEst, fileName_Table);
end

if(exportData2Excel)
    SelCoeff = selcCoeffNameAllItr_rep;
    Estimates = estimatesAllItr_rep;
    Method = methodsAllItr_rep;
    Valid = validAllItr_rep;
    TableOfSelecEst_rep = table(SelCoeff, Estimates, Method, Valid);
    fileName_Table = [dirNameScriptFile chosenSlash 'TableEstimates_5Sites_Rep' num2str(numRepComb) '_Set' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_numStrains' num2str(numStrainsInInitialPop) '_reg1' regStr1  '_reg2' regStr2 '.png' '.xlsx'];
    writetable(TableOfSelecEst_rep, fileName_Table);
end

%% plot histogram of number of estimated terms
% figure
% h1 = histogram(numSelcSitesEpiItr, [0.5:1:16.5])
% h1.FaceColor = color_scheme_npg(3,:);
% xlabel('Number of terms estimated')
% ylabel('Count')
% axis([0 16 0 numItr])
% % if(saveFigs == 1)
% %     figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_histNumTermsEst'];
% %     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 3.5])%6 4])%,[0 0 8 6.45])% ,[0 0 8 6])
% %     %set(gcf, 'renderer', 'painters');
% %     print([dirNameFigures figname], '-dpng','-r400')
% % end
% 
% figure
% h2 = histogram(numSelcSitesEpiItr_rep, [0.5:1:16.5])
% h2.FaceColor = color_scheme_npg(3,:);
% xlabel('Number of terms estimated')
% ylabel('Count')
% axis([0 16 0 numItr/5])
% if(saveFigs == 1)
%     figname = ['FiveSites_Rep' num2str(numRepComb) '_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_histNumTermsEst'];
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 3.5])%6 4])%,[0 0 8 6.45])% ,[0 0 8 6])
%     %set(gcf, 'renderer', 'painters');
%     print([dirNameFigures_rep figname], '-dpng','-r400')
% end




% one window bar plot
leftMargin = 0.25;
rightMargin = 0.02;
bottomMargin = 0.25;
topMargin = 0.1;
hgap = 0.0%0.09;
height2 = (1 - bottomMargin - topMargin);
height1 = height2;
width1 = (1 - leftMargin - rightMargin - hgap);
width2 = width1;


fig1 = figure
ha(1) = axes('Units','normalized', ...
                'Position',[leftMargin bottomMargin width1 height1], ...
                'XTickLabel','', ...
                'YTickLabel','');
            
            %
            
            
            
%LabelNumSiteSestSi = {'0.2', '0.4', '0.6', '0.8', '1'};
%LabelNumSiteSestSij = {'0.2', '0.4', '0.6', '0.8', '1'};
LabelNumSiteSestSi = {'0.2', '0.6', '1'};
LabelNumSiteSestSij = {'0.2', '0.6', '1'};

%h1 = histogram(numSelcSitesEpiItr, [0.5:1:16.5])
%h1 = histogram(numSelcSitesEpiItr_si, [0.5:1:16.5])
axes(ha(1))
h1 = histogram(numSelcSitesEpiItr_si_rep/5, [0:.2:1.2])
h1.FaceColor = color_scheme_npg(3,:);
set(gca, ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XGrid'       , 'off'      , ...
  'XColor'      , [.1 .1 .1], ...
  'YColor'      , [.1 .1 .1], ...
  'XTick', [0.3:0.4:1.2], ...
  'XTickLabel', LabelNumSiteSestSi, ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
  'LineWidth', 0.5)
ylabel('Count')
ylabh = get(gca,'ylabel'); 
set(ylabh,'Units','normalized');
set(ylabh,'position',get(ylabh,'position') + [-0.15 0 0]);

%xlabel('Number of accessible terms')
xlabel('Fraction of accessible terms')
xlabh = get(gca,'xlabel'); 
set(xlabh,'Units','normalized');
set(xlabh,'position',get(xlabh,'position') + [-0.15 -0.1 0]);
ylim([0 1000])
axis([0 1.2 0 numItr_rep])

hold on
%axes(ha(2))
h2 = histogram(numSelcSitesEpiItr_sij_rep/10, [0:.1:1.1])
h2.FaceColor = color_scheme_npg(5,:);
% set(gca, ...
%   'Box'         , 'on'     , ...
%   'TickDir'     , 'in'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'off'      , ...
%   'YMinorTick'  , 'off'      , ...
%   'YGrid'       , 'off'      , ...
%   'XGrid'       , 'off'      , ...
%   'XColor'      , [.1 .1 .1], ...
%   'YColor'      , [.1 .1 .1], ...
%   'XTick', [0.25:0.4:1.1], ...
%   'XTickLabel', LabelNumSiteSestSij, ... 
%   'YTickLabel', '', ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
%   'LineWidth', 0.5)
% axis([0 1.1 0 numItr])

if(saveFigs == 1)
    figname = ['FiveSites_Rep'  num2str(numRepComb) '_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_histNumTermsEst_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 3.5])%[0 0 6 4])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print([dirNameFigures_rep figname], '-dpng','-r400')
end


pause
close all
pause(2)

temp51 = strfind(fileName, 'itr');

fileNamePlotDataSave = ['PlotData_Rep'  num2str(numRepComb) '_' fileName(1:temp51+2) '1_' fileName(temp51+3:end)];
save([dirNameFigures_rep fileNamePlotDataSave])