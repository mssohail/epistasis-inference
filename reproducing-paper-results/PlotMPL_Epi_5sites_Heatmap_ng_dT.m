%%

% Last updated 03-Jan-2021
% code for 5-site simulation Figure 6

clc
close all
clear all

init_ColorNames_Settings;
% set(0,'DefaultAxesFontName','CMU Serif Roman')
% set(0,'DefaultTextFontName','CMU Serif Roman')
% set(0,'DefaultAxesFontSize',15)
% set(0,'DefaultTextFontSize',15)

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',7)
set(0,'DefaultTextFontSize',7)


exportData2Excel = 0;
saveFigs = 0;
saveFile = 1;
%Lin = 10;
thisSet = 1064701;%1062001%10603%1061003%10603%5611001%46%58%1991;%1990;%42%87%58%5%87%5;
fileNameContainingDirPath = 'dirNames.txt';
postFiltering = 1; % 0 : pre filtering (out put of AnalysisMPL_filt)
                   % 1 : post filtering
similarityThresh = 0.1; % threshold that sets allowed difference between two numbers to be still declared as same value, ie., 2 = 2.01
noiseThresh = 0.0;
getSysParam;
Tstart = 1;
Tused = 100;
ngAll = fliplr([20 50 100 200]);
dTAll = [5 10 20 25 50]%[5 10 20 30 50]% 100]; % [5 10 20 30 50 75]; % [5 10 20 30 50 100];
numStrainsInInitialPop = 20; % 5, 10, 30


dataFilesWith2Regs = 1;
if(dataFilesWith2Regs == 1)
     regStr1 = '1';
     regStr2 = '1';%'1';
end

numItr =1000%90%90%250;
%Tend = Tused + Tstart - 1;
Tend = Tused + Tstart;

actualT = T/1000;
posSelc = selVal/2/Nin;
negSelc = delSelVal/2/Nin;
lineCol = {'r.-', 'b.-', 'k.-', 'g.-'};
numTps = Tused/dT;
step = 0.002;
edges = [-10 -0.3:step:0.3 10] + step/2;

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
if(exist(dirNameFigures, 'dir') == 0)
    mkdir(dirNameFigures)        
end
selcCoeffNameColumnOneItr = ['s1  ';'s2  ';'s3  ';'s4  ';'s5  ';'ss12';'ss13';'ss14';'ss15';'ss23';'ss24';'ss25';'ss34';'ss35';'ss45'];
methodCol = [repmat('MPLE     ', 15,1) ; repmat('MPLE-NoMu', 15,1) ; repmat('MPL      ', 15,1) ; repmat('SL       ', 15,1)];



meanNormAbsErrHeatMap_MPL = zeros(length(ngAll), length(dTAll));
meanNormAbsErrHeatMap_cluster = zeros(length(ngAll), length(dTAll));

medianNormAbsErrHeatMap_MPL = zeros(length(ngAll), length(dTAll));
medianNormAbsErrHeatMap_cluster = zeros(length(ngAll), length(dTAll));


meanAUROC_ben_HeatMap_MPL = zeros(length(ngAll), length(dTAll));
meanAUROC_ben_HeatMap_Epi_si = zeros(length(ngAll), length(dTAll));
meanAUROC_ben_HeatMap_Epi_sij = zeros(length(ngAll), length(dTAll));
meanAUROC_ben_HeatMap_EpiWithR_si = zeros(length(ngAll), length(dTAll));
meanAUROC_ben_HeatMap_EpiWithR_sij = zeros(length(ngAll), length(dTAll));

meanAUROC_del_HeatMap_MPL = zeros(length(ngAll), length(dTAll));
meanAUROC_del_HeatMap_Epi_si = zeros(length(ngAll), length(dTAll));
meanAUROC_del_HeatMap_Epi_sij = zeros(length(ngAll), length(dTAll));
meanAUROC_del_HeatMap_EpiWithR_si = zeros(length(ngAll), length(dTAll));
meanAUROC_del_HeatMap_EpiWithR_sij = zeros(length(ngAll), length(dTAll));

meanAUROCSelc_ben_HeatMap_MPL = zeros(length(ngAll), length(dTAll));
meanAUROCSelc_ben_HeatMap_Epi_si = zeros(length(ngAll), length(dTAll));
meanAUROCSelc_ben_HeatMap_Epi_sij = zeros(length(ngAll), length(dTAll));
meanAUROCSelc_ben_HeatMap_EpiWithR_si = zeros(length(ngAll), length(dTAll));
meanAUROCSelc_ben_HeatMap_EpiWithR_sij = zeros(length(ngAll), length(dTAll));

meanAUROCSelc_del_HeatMap_MPL = zeros(length(ngAll), length(dTAll));
meanAUROCSelc_del_HeatMap_Epi_si = zeros(length(ngAll), length(dTAll));
meanAUROCSelc_del_HeatMap_Epi_sij = zeros(length(ngAll), length(dTAll));
meanAUROCSelc_del_HeatMap_EpiWithR_si = zeros(length(ngAll), length(dTAll));
meanAUROCSelc_del_HeatMap_EpiWithR_sij = zeros(length(ngAll), length(dTAll));

for nn = 1:length(ngAll)
    ng = ngAll(nn);
    for dd = 1:length(dTAll)
        dT = dTAll(dd);
        [ng dT]
selcCoeffNameAllItr = repmat(' ', numItr*60, 4);
estimatesAllItr = repmat(' ', numItr*60, 15);
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

aucLinkItr = -1*ones(numItr,2);
aucUnLinkItr = -1*ones(numItr,2);
aucLinkEpiItr = -1*ones(numItr,2);
aucLinkEpiWithRItr = -1*ones(numItr,2);
aucLinkNoMuEpiItr = -1*ones(numItr,2);
aucLinkEpiNoEItr = -1*ones(numItr,2);
aucLinkEpiOffTermsItr = -1*ones(numItr,2);
aucLinkNoMuEpiOffTermsItr = -1*ones(numItr,2);

aucLinkSelcItr = -1*ones(numItr,2);
aucUnLinkSelcItr = -1*ones(numItr,2);
aucLinkEpiSelcItr = -1*ones(numItr,2);
aucLinkEpiWithRSelcItr = -1*ones(numItr,2);
aucLinkNoMuEpiSelcItr = -1*ones(numItr,2);
aucLinkEpiNoESelcItr = -1*ones(numItr,2);
aucLinkEpiOffTermsSelcItr = -1*ones(numItr,2);
aucLinkNoMuEpiOffTermsSelcItr = -1*ones(numItr,2);

aucLinkEpi_only_si_Itr = -1*ones(numItr,2);
aucLinkEpi_only_sij_Itr = -1*ones(numItr,2);
aucLinkEpiWithR_only_si_Itr = -1*ones(numItr,2);
aucLinkEpiWithR_only_sij_Itr = -1*ones(numItr,2);
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

rhoItr = zeros(numItr,1);
numSelcSitesItr = zeros(1, numItr);
numSelcSitesEpiItr = zeros(1, numItr);
numSelcSitesEpiItr_si = zeros(1, numItr);
numSelcSitesEpiItr_sij = zeros(1, numItr);

allEstEpiItr = zeros(numItr, 15);
allEstEpiWithRItr = zeros(numItr, 15);
allEstNoMuEpiItr = zeros(numItr, 15);
allEstEpiNoEItr = zeros(numItr, 15);
allEstItr = zeros(numItr, 5);

allSelcSites = zeros(numItr, 5);
allSelcSitesEpi = zeros(numItr, 15);
for itr = 1:numItr%[1:50 52:100]%1:numItr
    
    selcSites = [];
    
%     if(dataFilesWith2Regs == 1)
%         if(Tstart == 1)
%             fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
%         else        
%             fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
%         end        
%     else
%         if(Tstart == 1)
%             fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
%         else        
%             fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
%         end
%     end
%     load([dirNameAnalysis fileName], 'L', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
%              'sigmaEstOutLink', 'sigmaEstOutLinkEpi', 'sigmaEstOutLinkNoMuEpi',...
%              'sigmaEstOutUnLink', 'sigmaEstOutLinkEpiNoE', ...
%              'q', 'q11', 'priorConst', 'timeWholeCodeMPL',...'perSiteSelctionSelc', 'selcSites', ... 'sigmaEstOutLinkSelc', 'sigmaEstOutLinkNoMuSelc','sigmaEstOutUnLinkSelc', 'sigmaEstOutUnLinkNoMuSelc', ...
%              'deltaLogLikeli', 'deltaLogLikeli50', 'deltaLogLikeli100',...
%              'deltaLogLikeliMPL', 'deltaLogLikeliSL', 'Hmat');
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
            else        
                fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
            end
        end
        load([dirNameAnalysis fileName], 'L', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
             'sigmaEstOutLink', 'sigmaEstOutLinkEpi', 'sigmaEstOutLinkEpiWithR', 'sigmaEstOutLinkNoMuEpi',...
             'sigmaEstOutUnLink', 'sigmaEstOutLinkEpiNoE', ...
             'q', 'q11', 'priorConst', 'timeWholeCodeMPL',...'perSiteSelctionSelc', 'selcSites', ... 'sigmaEstOutLinkSelc', 'sigmaEstOutLinkNoMuSelc','sigmaEstOutUnLinkSelc', 'sigmaEstOutUnLinkNoMuSelc', ...
             'deltaLogLikeli', 'deltaLogLikeli50', 'deltaLogLikeli100',...
             'deltaLogLikeliMPL', 'deltaLogLikeliSL', 'Hmat', ...
             'masterStrainList', 'currentNumOfCirStrainsInAllPop', 'masterSelectionFactor', 'freqStatsAllStrains');

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
             'sigmaEstOutUnLink', 'sigmaEstOutLinkEpiNoE', ...
             'q', 'q11', 'priorConst', 'timeWholeCodeMPL',...'perSiteSelctionSelc', 'selcSites', ... 'sigmaEstOutLinkSelc', 'sigmaEstOutLinkNoMuSelc','sigmaEstOutUnLinkSelc', 'sigmaEstOutUnLinkNoMuSelc', ...
             'deltaLogLikeli', 'deltaLogLikeli50', 'deltaLogLikeli100',...
             'deltaLogLikeliMPL', 'deltaLogLikeliSL', 'Hmat', ...
             'masterStrainList', 'currentNumOfCirStrainsInAllPop', 'masterSelectionFactor', 'freqStatsAllStrains');
             sigmaEstOutLinkEpiWithR = zeros(Lin*(Lin+1)/2, 1);
    end
%     deltaLL(itr) = deltaLogLikeli;
     
%     deltaLogLikeliMPLAll(itr) = deltaLogLikeliMPL;
%     deltaLogLikeliSLAll(itr) = deltaLogLikeliSL;

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
    %Identical_freqTraj;
    %[selcSitesEpi2, wellCondSites2, illCondSites2] = Identical_freqTraj(q, q11, Lin, similarityThresh, noiseThresh);
    [selcSitesEpi1, wellCondSites1, illCondSites1, cluster] = findIndColsOfHmat(Hmat*dT);
    %selcSitesEpi = selcSitesEpi1 & selcSitesEpi2;
    %wellCondSites = intersect(wellCondSites1, wellCondSites2);
    
    selcSitesEpi = selcSitesEpi1;
    wellCondSites = wellCondSites1;
    illCondSites = setdiff(1:15,wellCondSites);
    pairIDs = [1 2 3 4 5 12 13 14 15 23 24 25 34 35 45];
    selSitesEpiPairIDs{itr} = pairIDs(selcSitesEpi);
    
    if(postFiltering == 1)
        %selcSites = ~(abs(max(q) - min(q)) <= noiseThresh & (max(q) == 1 | min(q) == 0)); % ~ rejection criterion
        %selcSites = Identical_freqTraj_onlyMPL(q, Lin, similarityThresh, noiseThresh);
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
        sigmaEstOutLinkEpiWithRSelc = sigmaEstOutLinkEpiWithR(selcSitesEpi);
        sigmaEstOutLinkNoMuEpiSelc = sigmaEstOutLinkNoMuEpi(selcSitesEpi);
        sigmaEstOutLinkEpiNoESelc = sigmaEstOutLinkEpiNoE(selcSitesEpi);
        sigmaEstOutLinkSelc = sigmaEstOutLink(selcSites);
        sigmaEstOutUnLinkSelc = sigmaEstOutUnLink(selcSites);
        
%         if(noiseThresh == 0.05)
%             deltaLLSelc(itr) = deltaLogLikeli50;
%         elseif(noiseThresh == 0.1)
%             deltaLLSelc(itr) = deltaLogLikeli100;
%         else
%             disp('Error: no deltaLikelihood calculated for this noise thresh...')
%             pause
%         end
    end
    allSelcSites(itr,:) = selcSites;
    allSelcSitesEpi(itr,:) = selcSitesEpi;
    
    timeMPLItr(itr) = timeWholeCodeMPL(itr);
    % identify neutral, pos/neg selection

%     % needed for histograms
%     trueClassROC_AllItr = [trueClassROC_AllItr trueClassROC];
%     trueClassROCSelc_AllItr = [trueClassROCSelc_AllItr trueClassROCSelc];
%     trueClassROCEpi_AllItr = [trueClassROCEpi_AllItr trueClassROCEpi];
%     trueClassROCEpiSelc_AllItr = [trueClassROCEpiSelc_AllItr trueClassROCEpiSelc];
%     sigmaEstOutLinkEpi_AllItr = [sigmaEstOutLinkEpi_AllItr sigmaEstOutLinkEpi'];
%     sigmaEstOutLinkEpiSelc_AllItr = [sigmaEstOutLinkEpiSelc_AllItr sigmaEstOutLinkEpiSelc'];
%     
%     sigmaEstOutLinkNoMuEpi_AllItr = [sigmaEstOutLinkNoMuEpi_AllItr sigmaEstOutLinkNoMuEpi'];
%     sigmaEstOutLinkNoMuEpiSelc_AllItr = [sigmaEstOutLinkNoMuEpiSelc_AllItr sigmaEstOutLinkNoMuEpiSelc'];
%     
%     sigmaEstOutLinkEpiNoE_AllItr = [sigmaEstOutLinkEpiNoE_AllItr sigmaEstOutLinkEpiNoE'];
%     sigmaEstOutLinkEpiNoESelc_AllItr = [sigmaEstOutLinkEpiNoESelc_AllItr sigmaEstOutLinkEpiNoESelc'];
%     
%     sigmaEstOutLink_AllItr = [sigmaEstOutLink_AllItr sigmaEstOutLink'];
%     sigmaEstOutLinkSelc_AllItr = [sigmaEstOutLinkSelc_AllItr sigmaEstOutLinkSelc'];
%     
%     sigmaEstOutUnLink_AllItr = [sigmaEstOutUnLink_AllItr sigmaEstOutUnLink'];
%     sigmaEstOutUnLinkSelc_AllItr = [sigmaEstOutUnLinkSelc_AllItr sigmaEstOutUnLinkSelc'];
%     
%     perSiteSelction_AllItr = [perSiteSelction_AllItr perSiteSelction];
%     perSiteSelctionAllEpiTerms_AllItr = [perSiteSelctionAllEpiTerms_AllItr perSiteSelctionAllEpiTerms];

    
% % % % % % % %     nrmse_3class_Link_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutLink(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
% % % % % % % %     nrmse_3class_LinkNoMu_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutLinkNoMu(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
% % % % % % % %     nrmse_3class_UnLink_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutUnLink(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
% % % % % % % %     nrmse_3class_UnLinkNoMu_itrPos(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 4) - sigmaEstOutUnLinkNoMu(trueClassROC == 4)').^2)/sum(abs(perSiteSelction(trueClassROC == 4)).^2));
% % % % % % % % 
% % % % % % % %     nrmse_3class_LinkSelc_itrPos(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 4) - sigmaEstOutLinkSelc(trueClassROCSelc == 4)').^2)/sum(abs(perSiteSelctionSelc(trueClassROCSelc == 4)).^2));
% % % % % % % %     nrmse_3class_LinkNoMuSelc_itrPos(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 4) - sigmaEstOutLinkNoMuSelc(trueClassROCSelc == 4)').^2)/sum(abs(perSiteSelctionSelc(trueClassROCSelc == 4)).^2));
% % % % % % % %     nrmse_3class_UnLinkSelc_itrPos(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 4) - sigmaEstOutUnLinkSelc(trueClassROCSelc == 4)').^2)/sum(abs(perSiteSelctionSelc(trueClassROCSelc == 4)).^2));
% % % % % % % %     nrmse_3class_UnLinkNoMuSelc_itrPos(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 4) - sigmaEstOutUnLinkNoMuSelc(trueClassROCSelc == 4)').^2)/sum(abs(perSiteSelctionSelc(trueClassROCSelc == 4)).^2));
% % % % % % % % 
% % % % % % % %     nrmse_3class_Link_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutLink(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
% % % % % % % %     nrmse_3class_LinkNoMu_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutLinkNoMu(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
% % % % % % % %     nrmse_3class_UnLink_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutUnLink(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
% % % % % % % %     nrmse_3class_UnLinkNoMu_itrNeg(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 5) - sigmaEstOutUnLinkNoMu(trueClassROC == 5)').^2)/sum(abs(perSiteSelction(trueClassROC == 5)).^2));
% % % % % % % % 
% % % % % % % %     nrmse_3class_LinkSelc_itrNeg(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 5) - sigmaEstOutLinkSelc(trueClassROCSelc == 5)').^2)/sum(abs(perSiteSelctionSelc(trueClassROCSelc == 5)).^2));
% % % % % % % %     nrmse_3class_LinkNoMuSelc_itrNeg(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 5) - sigmaEstOutLinkNoMuSelc(trueClassROCSelc == 5)').^2)/sum(abs(perSiteSelctionSelc(trueClassROCSelc == 5)).^2));
% % % % % % % %     nrmse_3class_UnLinkSelc_itrNeg(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 5) - sigmaEstOutUnLinkSelc(trueClassROCSelc == 5)').^2)/sum(abs(perSiteSelctionSelc(trueClassROCSelc == 5)).^2));
% % % % % % % %     nrmse_3class_UnLinkNoMuSelc_itrNeg(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 5) - sigmaEstOutUnLinkNoMuSelc(trueClassROCSelc == 5)').^2)/sum(abs(perSiteSelctionSelc(trueClassROCSelc == 5)).^2));
% % % % % % % % 
% % % % % % % %     rmse_3class_Link_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutLink(trueClassROC == 3)').^2));
% % % % % % % %     rmse_3class_LinkNoMu_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutLinkNoMu(trueClassROC == 3)').^2));
% % % % % % % %     rmse_3class_UnLink_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutUnLink(trueClassROC == 3)').^2));
% % % % % % % %     rmse_3class_UnLinkNoMu_itrNeu(itr) = sqrt(sum(abs(perSiteSelction(trueClassROC == 3) - sigmaEstOutUnLinkNoMu(trueClassROC == 3)').^2));
% % % % % % % % 
% % % % % % % %     rmse_3class_LinkSelc_itrNeu(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 3) - sigmaEstOutLinkSelc(trueClassROCSelc == 3)').^2));
% % % % % % % %     rmse_3class_LinkNoMuSelc_itrNeu(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 3) - sigmaEstOutLinkNoMuSelc(trueClassROCSelc == 3)').^2));
% % % % % % % %     rmse_3class_UnLinkSelc_itrNeu(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 3) - sigmaEstOutUnLinkSelc(trueClassROCSelc == 3)').^2));
% % % % % % % %     rmse_3class_UnLinkNoMuSelc_itrNeu(itr) = sqrt(sum(abs(perSiteSelctionSelc(trueClassROCSelc == 3) - sigmaEstOutUnLinkNoMuSelc(trueClassROCSelc == 3)').^2));

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
    sitesEpi_only_si = logical([ones(1, 5), zeros(1, 10)]);
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
    sitesEpi_only_sij = logical([zeros(1, 5), ones(1, 10)]);
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
    % Calculate ABSOLUTE ERROR metrics
    if(sum(selcSitesEpi_only_si) > 0)
        if(sum(abs(perSiteSelction_only_si_Selc)) < zeroThresh)
            absErr_LinkEpi_only_si_Selc_itr(itr) = sum(abs(perSiteSelction_only_si_Selc - sigmaEstOutLinkEpi_only_si_Selc'))/(sum(abs(perSiteSelction_only_si_Selc)) + naeReg);%/length(perSiteSelction_only_si_Selc);
        else
            absErr_LinkEpi_only_si_Selc_itr(itr) = sum(abs(perSiteSelction_only_si_Selc - sigmaEstOutLinkEpi_only_si_Selc'))/sum(abs(perSiteSelction_only_si_Selc));%/length(perSiteSelction_only_si_Selc);
        end
    end
    if(sum(abs(perSiteSelction_only_si)) < zeroThresh)
        absErr_LinkEpi_only_si_itr(itr) = sum(abs(perSiteSelction_only_si - sigmaEstOutLinkEpi_only_si'))/(sum(abs(perSiteSelction_only_si)) + naeReg);%/length(perSiteSelction_only_si);
    else
        absErr_LinkEpi_only_si_itr(itr) = sum(abs(perSiteSelction_only_si - sigmaEstOutLinkEpi_only_si'))/sum(abs(perSiteSelction_only_si));%/length(perSiteSelction_only_si);
    end
    
    if(sum(selcSitesEpi_only_sij) > 0)
        if(sum(abs(perSiteSelction_only_sij_Selc)) < zeroThresh)
            absErr_LinkEpi_only_sij_Selc_itr(itr) = sum(abs(perSiteSelction_only_sij_Selc - sigmaEstOutLinkEpi_only_sij_Selc'))/(sum(abs(perSiteSelction_only_sij_Selc)) + naeReg);%/length(perSiteSelction_only_sij_Selc);
        else
            absErr_LinkEpi_only_sij_Selc_itr(itr) = sum(abs(perSiteSelction_only_sij_Selc - sigmaEstOutLinkEpi_only_sij_Selc'))/sum(abs(perSiteSelction_only_sij_Selc));%/length(perSiteSelction_only_sij_Selc);
        end
    end
    if(sum(abs(perSiteSelction_only_sij)) < zeroThresh)
        absErr_LinkEpi_only_sij_itr(itr) = sum(abs(perSiteSelction_only_sij - sigmaEstOutLinkEpi_only_sij'))/(sum(abs(perSiteSelction_only_sij)) + naeReg);%/length(perSiteSelction_only_sij);
    else
        absErr_LinkEpi_only_sij_itr(itr) = sum(abs(perSiteSelction_only_sij - sigmaEstOutLinkEpi_only_sij'))/sum(abs(perSiteSelction_only_sij));%/length(perSiteSelction_only_sij);
    end
    
    if(sum(selcSitesEpi) > 0)
        if(abs(perSiteSelctionAllEpiTerms(selcSitesEpi)) < zeroThresh)
            absErr_LinkEpi_Selc_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi) - sigmaEstOutLinkEpiSelc'))/(sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi))) + naeReg);%/length(perSiteSelctionAllEpiTerms(selcSitesEpi));
        else
            absErr_LinkEpi_Selc_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi) - sigmaEstOutLinkEpiSelc'))/sum(abs(perSiteSelctionAllEpiTerms(selcSitesEpi)));%/length(perSiteSelctionAllEpiTerms(selcSitesEpi));
        end
    end
    if(sum(abs(perSiteSelctionAllEpiTerms)) < zeroThresh)
        absErr_LinkEpi_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms - sigmaEstOutLinkEpi'))/(sum(abs(perSiteSelctionAllEpiTerms)) + naeReg);%/length(perSiteSelctionAllEpiTerms);
    else
        absErr_LinkEpi_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms - sigmaEstOutLinkEpi'))/sum(abs(perSiteSelctionAllEpiTerms));%/length(perSiteSelctionAllEpiTerms);
    end
    
    if(sum(selcSites) > 0)
        if(sum(abs(perSiteSelction(selcSites))) < zeroThresh)
            absErr_LinkSelc_itr(itr) = sum(abs(perSiteSelction(selcSites) - sigmaEstOutLinkSelc'))/(sum(abs(perSiteSelction(selcSites))) + naeReg);%/length(perSiteSelction(selcSites));
            absErr_UnLinkSelc_itr(itr) = sum(abs(perSiteSelction(selcSites) - sigmaEstOutUnLinkSelc'))/(sum(abs(perSiteSelction(selcSites))) + naeReg);%/length(perSiteSelction(selcSites));
        else
            absErr_LinkSelc_itr(itr) = sum(abs(perSiteSelction(selcSites) - sigmaEstOutLinkSelc'))/sum(abs(perSiteSelction(selcSites)));%/length(perSiteSelction(selcSites));
            absErr_UnLinkSelc_itr(itr) = sum(abs(perSiteSelction(selcSites) - sigmaEstOutUnLinkSelc'))/sum(abs(perSiteSelction(selcSites)));%/length(perSiteSelction(selcSites));
        end
    end
    if(sum(abs(perSiteSelction)) < zeroThresh)
        absErr_Link_itr(itr) = sum(abs(perSiteSelction - sigmaEstOutLink'))/(sum(abs(perSiteSelction)) + naeReg);%/length(perSiteSelction);
        absErr_UnLink_itr(itr) = sum(abs(perSiteSelction - sigmaEstOutUnLink'))/(sum(abs(perSiteSelction)) + naeReg);%/length(perSiteSelction);
    else
        absErr_Link_itr(itr) = sum(abs(perSiteSelction - sigmaEstOutLink'))/sum(abs(perSiteSelction));%/length(perSiteSelction);
        absErr_UnLink_itr(itr) = sum(abs(perSiteSelction - sigmaEstOutUnLink'))/sum(abs(perSiteSelction));%/length(perSiteSelction);
    end
    
    perSiteSelctionAllEpiTerms_cluster = perSiteSelctionAllEpiTerms(wellCondSites);
    sigmaEstOutLinkEpi_cluster = sigmaEstOutLinkEpi(wellCondSites)';
    perSiteSelctionAllEpiTerms_clusterComb = perSiteSelctionAllEpiTerms_cluster;
    sigmaEstOutLinkEpi_clusterComb = sigmaEstOutLinkEpi_cluster;
    indTemp201 = [];
    num_clusterWellOnly = length(sigmaEstOutLinkEpi_cluster);
    num_clusterCombWellOnly = length(sigmaEstOutLinkEpi_clusterComb);
    for wt = 1:length(cluster)
        indTemp101 = cluster{wt};
        indTemp201 = [indTemp201 indTemp101];
        perSiteSelctionAllEpiTerms_cluster = [perSiteSelctionAllEpiTerms_cluster sum(perSiteSelctionAllEpiTerms(indTemp101))];
        sigmaEstOutLinkEpi_cluster = [sigmaEstOutLinkEpi_cluster sum(sigmaEstOutLinkEpi(indTemp101))];
    end
    clustInd = zeros(1, length(sigmaEstOutLinkEpi_cluster));
    temp65 = 1:num_clusterWellOnly;
    clustInd(temp65) = 1;
    clustInd = logical(clustInd);
    
    perSiteSelctionAllEpiTerms_clusterComb = [perSiteSelctionAllEpiTerms_clusterComb sum(perSiteSelctionAllEpiTerms(indTemp201))];
    sigmaEstOutLinkEpi_clusterComb = [sigmaEstOutLinkEpi_clusterComb sum(sigmaEstOutLinkEpi(indTemp201))];
    clustCombInd = zeros(1, length(sigmaEstOutLinkEpi_clusterComb));
    temp66 = 1:num_clusterCombWellOnly;
    clustCombInd(temp66) = 1;
    clustCombInd = logical(clustCombInd);
    
    
    perSiteSelctionAllEpiTerms_cluster = round(100000*perSiteSelctionAllEpiTerms_cluster)/100000;
    if(sum(abs(perSiteSelctionAllEpiTerms_cluster)) < zeroThresh)
        absErr_LinkEpi_cluster_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_cluster - sigmaEstOutLinkEpi_cluster))/sum(abs(perSiteSelctionAllEpiTerms_cluster) + naeReg);%/length(perSiteSelctionAllEpiTerms_cluster);
    else
        absErr_LinkEpi_cluster_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_cluster - sigmaEstOutLinkEpi_cluster))/sum(abs(perSiteSelctionAllEpiTerms_cluster));%/length(perSiteSelctionAllEpiTerms_cluster);
    end
    numTerms_LinkEpi_cluster_itr(itr) = length(perSiteSelctionAllEpiTerms_cluster);
    
    cluster_detail_cell{itr,1} = perSiteSelctionAllEpiTerms_cluster;
    cluster_detail_cell{itr,2} = sigmaEstOutLinkEpi_cluster;
    
    if(sum(clustInd) > 0)
        if(sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd))) < zeroThresh)
            absErr_LinkEpi_clusterOnlyWell_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd) - sigmaEstOutLinkEpi_cluster(clustInd)))/(sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd))) + naeReg);%/sum(clustInd);
        else
            absErr_LinkEpi_clusterOnlyWell_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd) - sigmaEstOutLinkEpi_cluster(clustInd)))/sum(abs(perSiteSelctionAllEpiTerms_cluster(clustInd)));%/sum(clustInd);
        end
    end
    
    if(sum(~clustInd) > 0)
        if(sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd))) < zeroThresh)
            absErr_LinkEpi_clusterOnlyAmb_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd) - sigmaEstOutLinkEpi_cluster(~clustInd)))/(sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd))) + naeReg);%/sum(~clustInd);
        else
            absErr_LinkEpi_clusterOnlyAmb_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd) - sigmaEstOutLinkEpi_cluster(~clustInd)))/sum(abs(perSiteSelctionAllEpiTerms_cluster(~clustInd)));%/sum(~clustInd);
        end
    end
    
    if(isinf(absErr_LinkEpi_cluster_itr(itr)))
        pause
    end
    
    if(absErr_LinkEpi_cluster_itr(itr) > 5)
        %pause
    end
    
    % cluster     : well cond sites unique and ill cond sites divided into 
    %               amb clusters, each cluster as 1 selCoeff
    % clusterComb : all amb cluster bunched in 1 selCoeff
    
    
    perSiteSelctionAllEpiTerms_clusterComb = round(100000*perSiteSelctionAllEpiTerms_clusterComb)/100000;
    if(sum(abs(perSiteSelctionAllEpiTerms_clusterComb)) < zeroThresh)
        absErr_LinkEpi_clusterComb_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_clusterComb - sigmaEstOutLinkEpi_clusterComb))/sum(abs(perSiteSelctionAllEpiTerms_clusterComb) + naeReg);%/length(perSiteSelctionAllEpiTerms_clusterComb);
    else
        absErr_LinkEpi_clusterComb_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_clusterComb - sigmaEstOutLinkEpi_clusterComb))/sum(abs(perSiteSelctionAllEpiTerms_clusterComb));%/length(perSiteSelctionAllEpiTerms_clusterComb);
    end
    numTerms_LinkEpi_clusterComb_itr(itr) = length(perSiteSelctionAllEpiTerms_clusterComb);
    
    if(sum(clustCombInd) > 0)
        if(sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd))) < zeroThresh)
            absErr_LinkEpi_clusterCombOnlyWell_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd) - sigmaEstOutLinkEpi_clusterComb(clustCombInd)))/(sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd))) + naeReg);%/sum(clustCombInd);
        else
            absErr_LinkEpi_clusterCombOnlyWell_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd) - sigmaEstOutLinkEpi_clusterComb(clustCombInd)))/sum(abs(perSiteSelctionAllEpiTerms_clusterComb(clustCombInd)));%/sum(clustCombInd);
        end
    end
    
    if(sum(~clustCombInd) > 0)
        if(sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd))) < zeroThresh)
            absErr_LinkEpi_clusterCombOnlyAmb_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd) - sigmaEstOutLinkEpi_clusterComb(~clustCombInd)))/(sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd))) + naeReg);%/sum(~clustCombInd);
        else
            absErr_LinkEpi_clusterCombOnlyAmb_itr(itr) = sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd) - sigmaEstOutLinkEpi_clusterComb(~clustCombInd)))/sum(abs(perSiteSelctionAllEpiTerms_clusterComb(~clustCombInd)));%/sum(~clustCombInd);
        end
    end
    
    
    if(absErr_LinkEpi_clusterComb_itr(itr) > 40)
%         pause
    end
%     selcSitesEpi_only_sij = selcSitesEpi & logical([zeros(1, 5), ones(1, 10)]);
%     perSiteSelction_only_sij_Selc = perSiteSelctionAllEpiTerms(selcSitesEpi_only_sij);
%     sigmaEstOutLinkEpi_only_sij_Selc = sigmaEstOutLinkEpi(selcSitesEpi_only_sij);
% 
%     posOnlyItrTemp_only_sij = posOnlyItrTemp_AllEpiTerms & logical([zeros(1, 5), ones(1, 10)]);
%     negOnlyItrTemp_only_sij = negOnlyItrTemp_AllEpiTerms & logical([zeros(1, 5), ones(1, 10)]);
%     posOnlySelcItrTemp_only_sij = perSiteSelction_only_sij_Selc > neutralLimit;
%     negOnlySelcItrTemp_only_sij = perSiteSelction_only_sij_Selc < -neutralLimit;
%     [aucLinkEpiItrTemp_only_sij, aucLinkEpiSelcItrTemp_only_sij] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelction_only_sij_Selc, sigmaEstOutLinkEpi, sigmaEstOutLinkEpi_only_sij_Selc, posOnlyItrTemp_only_sij, negOnlyItrTemp_only_sij, posOnlySelcItrTemp_only_sij, negOnlySelcItrTemp_only_sij);
%     %[aucLinkNoMuEpiItrTemp, aucLinkNoMuEpiSelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkNoMuEpi, sigmaEstOutLinkNoMuEpiSelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);
%     %[aucLinkEpiNoEItrTemp, aucLinkEpiNoESelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkEpiNoE, sigmaEstOutLinkEpiNoESelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);

    
    % Save all variables for this itr
    %--------------------------------------

    % Save AUROC unfiltered
    aucLinkEpiItr(itr,:) = aucLinkEpiItrTemp;
    aucLinkEpiWithRItr(itr,:) = aucLinkEpiWithRItrTemp;
    aucLinkNoMuEpiItr(itr,:) = aucLinkNoMuEpiItrTemp;
    aucLinkEpiNoEItr(itr,:) = aucLinkEpiNoEItrTemp;
    aucLinkItr(itr,:) = aucLinkItrTemp;
    aucUnLinkItr(itr,:) = aucUnLinkItrTemp;

    % Save AUROC filtered
    aucLinkEpiSelcItr(itr,:) = aucLinkEpiSelcItrTemp;
    aucLinkEpiWithRSelcItr(itr,:) = aucLinkEpiWithRSelcItrTemp;
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
    allEstEpiWithRItr(itr,:) = sigmaEstOutLinkEpiWithR';
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

       meanNormAbsErrHeatMap_MPL(nn,dd) = meanAllNormAbsErrs_Selc(6);
       meanNormAbsErrHeatMap_cluster(nn,dd) = meanAllNormAbsErrs(4);
       
       medianNormAbsErrHeatMap_MPL(nn,dd) = medianAllNormAbsErrs_Selc(6);
       medianNormAbsErrHeatMap_cluster(nn,dd) = medianAllNormAbsErrs(4);
       
            % disp('-------------------------------')
            % disp(' mean AUROC   /     AUROCSelc')
            % allAuc = [mean(aucLinkEpiItr) mean(aucLinkEpiSelcItr);
            % mean(aucLinkNoMuEpiItr) mean(aucLinkNoMuEpiSelcItr);
            % %mean(aucLinkEpiNoEItr) mean(aucLinkEpiNoESelcItr);
            % mean(aucLinkItr) mean(aucLinkSelcItr);
            % mean(aucUnLinkItr) mean(aucUnLinkSelcItr)]


            
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
            mean(aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,2) ~= -1,2))]
            disp('-------------------------------')       

            
       meanAUROC_ben_HeatMap_MPL(nn,dd) = allAuc(3,1);
       meanAUROC_ben_HeatMap_Epi_si(nn,dd) = allAuc(5,1);
       meanAUROC_ben_HeatMap_Epi_sij(nn,dd) = allAuc(6,1);
       meanAUROC_ben_HeatMap_EpiWithR_si(nn,dd) = allAuc(7,1);
       meanAUROC_ben_HeatMap_EpiWithR_sij(nn,dd) = allAuc(8,1);
       meanAUROCSelc_ben_HeatMap_MPL(nn,dd) = allAuc(3,3);
       meanAUROCSelc_ben_HeatMap_Epi_si(nn,dd) = allAuc(5,3);
       meanAUROCSelc_ben_HeatMap_Epi_sij(nn,dd) = allAuc(6,3);
       meanAUROCSelc_ben_HeatMap_EpiWithR_si(nn,dd) = allAuc(7,3);
       meanAUROCSelc_ben_HeatMap_EpiWithR_sij(nn,dd) = allAuc(8,3);
       
       meanAUROC_del_HeatMap_MPL(nn,dd) = allAuc(3,2);
       meanAUROC_del_HeatMap_Epi_si(nn,dd) = allAuc(5,2);
       meanAUROC_del_HeatMap_Epi_sij(nn,dd) = allAuc(6,2);
       meanAUROC_del_HeatMap_EpiWithR_si(nn,dd) = allAuc(7,2);
       meanAUROC_del_HeatMap_EpiWithR_sij(nn,dd) = allAuc(8,2);
       meanAUROCSelc_del_HeatMap_MPL(nn,dd) = allAuc(3,4);
       meanAUROCSelc_del_HeatMap_Epi_si(nn,dd) = allAuc(5,4);
       meanAUROCSelc_del_HeatMap_Epi_sij(nn,dd) = allAuc(6,4);
       meanAUROCSelc_del_HeatMap_EpiWithR_si(nn,dd) = allAuc(7,4);
       meanAUROCSelc_del_HeatMap_EpiWithR_sij(nn,dd) = allAuc(8,4);
    end
end

%%
if(saveFile == 1)
    fileNameSave = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_HeatmapNgdT_ng' num2str(ngAll(end)) '_' num2str(ngAll(1)) '_dT' num2str(dTAll(1)) '_' num2str(dTAll(end)) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '_numItr' num2str(numItr) '.mat'];
    save([dirNameFigures fileNameSave])
end