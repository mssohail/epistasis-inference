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
thisSetAll = 1064701%[1066201:10:1066291 1066401:10:1066491 1066601:10:1066691 1066801:10:1066891]
for setInd = 1:length(thisSetAll)
    thisSet = thisSetAll(setInd);%1068242%[1068141 1068144 1068241 1068244 1068143 1068142 1068243 1068242]
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
    numStrainsInInitialPop = 20%5; % 5, 10, 30

    dataFilesWith2Regs = 1;
    if(dataFilesWith2Regs == 1)
         regStr1 = '1';
         regStr2 = '1';
    end

    thisItrStart = 1;
    thisItrEnd = 5;%5;
    lastItr = 1000;%1000;
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
    sigmaEstOutLinkEpiWithR_AllItr = [];
    sigmaEstOutLinkEpiWithRSelc_AllItr = [];
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
    aucLinkEpiWithRItr = -1*ones(numItr,2);
    aucLinkNoMuEpiItr = -1*ones(numItr,2);
    aucLinkEpiNoEItr = -1*ones(numItr,2);
    aucLinkItr = -1*ones(numItr,2);
    aucUnLinkItr = -1*ones(numItr,2);
    aucLinkEpiSelcItr = -1*ones(numItr,2);
    aucLinkEpiWithRSelcItr = -1*ones(numItr,2);
    aucLinkNoMuEpiSelcItr = -1*ones(numItr,2);
    aucLinkEpiNoESelcItr = -1*ones(numItr,2);
    aucLinkSelcItr = -1*ones(numItr,2);
    aucUnLinkSelcItr = -1*ones(numItr,2);

    numSelcSitesItr = zeros(1, numItr);
    numSelcSitesEpiItr = zeros(1, numItr);
    numSelcSitesEpiItr_si = zeros(1, numItr);
    numSelcSitesEpiItr_sij = zeros(1, numItr);
    numSelcSitesEpiWithRItr_si = zeros(1, numItr);
    numSelcSitesEpiWithRItr_sij = zeros(1, numItr);

    allEstEpiItr = zeros(numItr, numParam);
    allEstEpiWithRItr = zeros(numItr, numParam);
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
    aucLinkEpiWithRItr_rep = -1*ones(numItr_rep,2);
    aucLinkNoMuEpiItr_rep = -1*ones(numItr_rep,2);
    aucLinkEpiNoEItr_rep = -1*ones(numItr_rep,2);
    aucLinkItr_rep = -1*ones(numItr_rep,2);
    aucUnLinkItr_rep = -1*ones(numItr_rep,2);
    aucLinkEpiSelcItr_rep = -1*ones(numItr_rep,2);
    aucLinkEpiWithRSelcItr_rep = -1*ones(numItr_rep,2);
    aucLinkNoMuEpiSelcItr_rep = -1*ones(numItr_rep,2);
    aucLinkEpiNoESelcItr_rep = -1*ones(numItr_rep,2);
    aucLinkSelcItr_rep = -1*ones(numItr_rep,2);
    aucUnLinkSelcItr_rep = -1*ones(numItr_rep,2);

    numSelcSitesItr_rep = zeros(1, numItr_rep);
    numSelcSitesEpiItr_rep = zeros(1, numItr_rep);
    numSelcSitesEpiItr_si_rep = zeros(1, numItr_rep);
    numSelcSitesEpiItr_sij_rep = zeros(1, numItr_rep);
    numSelcSitesEpiWithRItr_si_rep = zeros(1, numItr_rep);
    numSelcSitesEpiWithRItr_sij_rep = zeros(1, numItr_rep);

    allEstEpiItr_rep = zeros(numItr_rep, numParam);
    allEstEpiWithRItr_rep = zeros(numItr_rep, numParam);
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
            fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
        else        
            fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
        end
        load([dirNameAnalysis fileName], 'Lin', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
                 'sigmaEstOutLink', 'sigmaEstOutLinkEpi', 'sigmaEstOutLinkNoMuEpi',...
                 'sigmaEstOutUnLink', 'sigmaEstOutLinkEpiNoE', 'sigmaEstOutLinkEpiWithR', ...
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
            sigmaEstOutLinkEpiWithRSelc = sigmaEstOutLinkEpiWithR(selcSitesEpi);
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
        aucLinkEpiWithR_only_si_Itr(itr,:) = aucLinkEpiWithRItrTemp_only_si;
        aucLinkEpiWithR_only_sij_Itr(itr,:) = aucLinkEpiWithRItrTemp_only_sij;

        aucLinkEpi_only_si_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_si;
        aucLinkEpi_only_sij_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_sij;
        aucLinkEpiWithR_only_si_SelcItr(itr,:) = aucLinkEpiWithRSelcItrTemp_only_si;
        aucLinkEpiWithR_only_sij_SelcItr(itr,:) = aucLinkEpiWithRSelcItrTemp_only_sij;


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

        %jj = jj + 1;
    %-------------------------------------------------------------------------------------------    
        % 1.2 Load data (combining replicates)
        if(rem(itr, thisItrEnd) == 0)    

            if(Tstart == 1)
                fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff'  '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str((xCount-1)*numRepComb +thisItrStart) '_' num2str(xCount*numRepComb) '.mat'];
            else        
                fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str((xCount-1)*numRepComb +thisItrStart) '_' num2str(xCount*numRepComb) '.mat'];
            end

            load([dirNameAnalysis fileName], 'Lin', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
                 'sigmaEstOutLinkRep', 'sigmaEstOutLinkEpiRep', 'sigmaEstOutLinkEpiWithRRep',...
                 'sigmaEstOutLinkNoMuEpiRep',...
                 'sigmaEstOutUnLinkRep', 'sigmaEstOutLinkEpiNoERep', 'HmatDTAll');


            sigmaEstOutLink_rep = sigmaEstOutLinkRep;
            sigmaEstOutLinkEpi_rep = sigmaEstOutLinkEpiRep;
            sigmaEstOutLinkEpiWithR_rep = sigmaEstOutLinkEpiWithRRep;
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
            sigmaEstOutLinkEpiWithRSelc_rep = sigmaEstOutLinkEpiWithR(selcSitesEpi_rep);
            sigmaEstOutLinkNoMuEpiSelc_rep = sigmaEstOutLinkNoMuEpi(selcSitesEpi_rep);
            sigmaEstOutLinkEpiNoESelc_rep = sigmaEstOutLinkEpiNoE(selcSitesEpi_rep);
            sigmaEstOutLinkSelc_rep = sigmaEstOutLink(selcSites_rep);
            sigmaEstOutUnLinkSelc_rep = sigmaEstOutUnLink(selcSites_rep);




            [aucLinkItrTemp_rep, aucLinkSelcItrTemp_rep, aucUnLinkItrTemp_rep, aucUnLinkSelcItrTemp_rep, ...
            aucLinkEpiItrTemp_rep, aucLinkEpiSelcItrTemp_rep, aucLinkNoMuEpiItrTemp_rep, aucLinkNoMuEpiSelcItrTemp_rep, ...    
            aucLinkEpiNoEItrTemp_rep, aucLinkEpiNoESelcItrTemp_rep, aucLinkEpiItrTemp_only_si_rep, aucLinkEpiSelcItrTemp_only_si_rep, ...
            aucLinkEpiWithRItrTemp_only_si_rep, aucLinkEpiWithRSelcItrTemp_only_si_rep, aucLinkEpiItrTemp_only_sij_rep, aucLinkEpiSelcItrTemp_only_sij_rep, ...
            aucLinkEpiWithRItrTemp_only_sij_rep, aucLinkEpiWithRSelcItrTemp_only_sij_rep] = calcEpiWithRRef_AUROC_new(perSiteSelction, perSiteSelctionSelc_rep, ...
                                perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc_rep, selcSitesEpi_rep, ...
                                sigmaEstOutLink_rep, sigmaEstOutLinkSelc_rep, sigmaEstOutUnLink_rep, sigmaEstOutUnLinkSelc_rep, ...
                                sigmaEstOutLinkEpi_rep, sigmaEstOutLinkEpiSelc_rep, sigmaEstOutLinkNoMuEpi_rep, sigmaEstOutLinkNoMuEpiSelc_rep, ...
                                sigmaEstOutLinkEpiNoE_rep, sigmaEstOutLinkEpiNoESelc_rep, neutralLimit, ...
                                sigmaEstOutLinkEpiWithR_rep, sigmaEstOutLinkEpiWithRSelc_rep);

            [aucLinkItrTemp, aucLinkSelcItrTemp, aucUnLinkItrTemp, aucUnLinkSelcItrTemp, ...
            aucLinkEpiItrTemp, aucLinkEpiSelcItrTemp, aucLinkNoMuEpiItrTemp, aucLinkNoMuEpiSelcItrTemp, ...    
            aucLinkEpiNoEItrTemp, aucLinkEpiNoESelcItrTemp, aucLinkEpiItrTemp_only_si, aucLinkEpiSelcItrTemp_only_si, ...
            aucLinkEpiWithRItrTemp_only_si, aucLinkEpiWithRSelcItrTemp_only_si, aucLinkEpiItrTemp_only_sij, aucLinkEpiSelcItrTemp_only_sij, ...
            aucLinkEpiWithRItrTemp_only_sij, aucLinkEpiWithRSelcItrTemp_only_sij] =  calcEpiWithRRef_AUROC_new(perSiteSelction, perSiteSelctionSelc, ...
                                        perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, selcSitesEpi, ...
                                        sigmaEstOutLink, sigmaEstOutLinkSelc, sigmaEstOutUnLink, sigmaEstOutUnLinkSelc, ...
                                        sigmaEstOutLinkEpi, sigmaEstOutLinkEpiSelc, sigmaEstOutLinkNoMuEpi, sigmaEstOutLinkNoMuEpiSelc, ...
                                        sigmaEstOutLinkEpiNoE, sigmaEstOutLinkEpiNoESelc, neutralLimit, ...
                                        sigmaEstOutLinkEpiWithR, sigmaEstOutLinkEpiWithRSelc);


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
            aucLinkEpiWithR_only_si_Itr_rep(xCount,:) = aucLinkEpiWithRItrTemp_only_si_rep;
            aucLinkEpiWithR_only_sij_Itr_rep(xCount,:) = aucLinkEpiWithRItrTemp_only_sij_rep;

            aucLinkEpi_only_si_SelcItr_rep(xCount,:) = aucLinkEpiSelcItrTemp_only_si_rep;
            aucLinkEpi_only_sij_SelcItr_rep(xCount,:) = aucLinkEpiSelcItrTemp_only_sij_rep;
            aucLinkEpiWithR_only_si_SelcItr_rep(xCount,:) = aucLinkEpiWithRSelcItrTemp_only_si_rep;
            aucLinkEpiWithR_only_sij_SelcItr_rep(xCount,:) = aucLinkEpiWithRSelcItrTemp_only_sij_rep;

            condHmatWellItr_rep(itr) = cond(HmatDTAll(wellCondSites_rep,wellCondSites_rep));

            allEstEpiItr_rep(xCount,:) = sigmaEstOutLinkEpi_rep';
            allEstEpiWithRItr_rep(xCount,:) = sigmaEstOutLinkEpiWithR_rep';
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

    % %%
    % disp('-------------------------------')
    % disp(' mean AUROC   /     AUROCSelc')
    % allAuc = [mean(aucLinkEpiItr(aucLinkEpiItr(:,1) ~= -1, 1)) mean(aucLinkEpiItr(aucLinkEpiItr(:,2) ~= -1, 2)) mean(aucLinkEpiSelcItr(aucLinkEpiSelcItr(:,1) ~= -1,1)) mean(aucLinkEpiSelcItr(aucLinkEpiSelcItr(:,2) ~= -1,2));
    % mean(aucLinkNoMuEpiItr(aucLinkNoMuEpiItr(:,1) ~= -1, 1)) mean(aucLinkNoMuEpiItr(aucLinkNoMuEpiItr(:,2) ~= -1, 2)) mean(aucLinkNoMuEpiSelcItr(aucLinkNoMuEpiSelcItr(:,1) ~= -1,1)) mean(aucLinkNoMuEpiSelcItr(aucLinkNoMuEpiSelcItr(:,2) ~= -1,2));
    % %mean(aucLinkEpiNoEItr) mean(aucLinkEpiNoESelcItr);
    % mean(aucLinkItr(aucLinkItr(:,1) ~= -1, 1)) mean(aucLinkItr(aucLinkItr(:,2) ~= -1, 2)) mean(aucLinkSelcItr(aucLinkSelcItr(:,1) ~= -1,1)) mean(aucLinkSelcItr(aucLinkSelcItr(:,2) ~= -1,2));
    % mean(aucUnLinkItr(aucUnLinkItr(:,1) ~= -1, 1)) mean(aucUnLinkItr(aucUnLinkItr(:,2) ~= -1, 2)) mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,1) ~= -1,1)) mean(aucUnLinkSelcItr(aucUnLinkSelcItr(:,2) ~= -1,2));
    % mean(aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_Itr(aucLinkEpi_only_si_Itr(:,2) ~= -1,2)) mean(aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpi_only_si_SelcItr(aucLinkEpi_only_si_SelcItr(:,2) ~= -1,2));
    % mean(aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,2) ~= -1,2)) mean(aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,2) ~= -1,2));
    % mean(aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,2) ~= -1,2));
    % mean(aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,2) ~= -1,2))]
    % disp('-------------------------------')

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
    mean(aucLinkEpi_only_sij_Itr_rep(aucLinkEpi_only_sij_Itr_rep(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_Itr_rep(aucLinkEpi_only_sij_Itr_rep(:,2) ~= -1,2)) mean(aucLinkEpi_only_sij_SelcItr_rep(aucLinkEpi_only_sij_SelcItr_rep(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_SelcItr_rep(aucLinkEpi_only_sij_SelcItr_rep(:,2) ~= -1,2));
    mean(aucLinkEpiWithR_only_si_Itr_rep(aucLinkEpiWithR_only_si_Itr_rep(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_Itr_rep(aucLinkEpiWithR_only_si_Itr_rep(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_si_SelcItr_rep(aucLinkEpiWithR_only_si_SelcItr_rep(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_SelcItr_rep(aucLinkEpiWithR_only_si_SelcItr_rep(:,2) ~= -1,2));
    mean(aucLinkEpiWithR_only_sij_Itr_rep(aucLinkEpiWithR_only_sij_Itr_rep(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_Itr_rep(aucLinkEpiWithR_only_sij_Itr_rep(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_sij_SelcItr_rep(aucLinkEpiWithR_only_sij_SelcItr_rep(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_SelcItr_rep(aucLinkEpiWithR_only_sij_SelcItr_rep(:,2) ~= -1,2))]
    disp('-------------------------------')




    %%

    temp51 = strfind(fileName, 'itr');

    fileNamePlotDataSave = ['PlotData_Rep'  num2str(numRepComb) '_' fileName(1:temp51+2) '1_' fileName(temp51+3:end)];
    save([dirNameFigures_rep fileNamePlotDataSave])
end