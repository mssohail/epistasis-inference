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

allSets = [1065401 1066201:10:1066291 1066401:10:1066491 1066601:10:1066691 1066801:10:1066891 1064701];

allNumItr = [1000 100*ones(1, 40) 1000];


for ii = 1:length(allSets)

    thisSet = allSets(ii);%1066891;%1065301;%1044001;%1064701;%1064501;%1064001;
    numItr = allNumItr(ii)%90%90%90;

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
    numStrainsInInitialPop = 20%10;%5; % 5, 10, 30
    useStrongReg = 0;

    dataFilesWith2Regs = 1;
    if(dataFilesWith2Regs == 1)
        regStr1 = '1';
        regStr2 = '1';
    end


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
    aucLinkEpi_only_sij_Itr = -1*ones(numItr,2);
    aucLinkEpi_only_si_SelcItr = -1*ones(numItr,2);
    aucLinkEpi_only_sij_SelcItr = -1*ones(numItr,2);
    aucLinkEpiWithR_only_si_SelcItr = -1*ones(numItr,2);

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

        aucLinkEpi_only_si_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_si;
        aucLinkEpi_only_sij_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_sij;
        aucLinkEpiWithR_only_si_SelcItr(itr,:) = aucLinkEpiWithRSelcItrTemp_only_si;

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
    mean(aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,2) ~= -1,2));]
    disp('-------------------------------')

    %%
    % bar plot of MPLE and MPL
    %----------------------------
    for k = 1:2
        fig1 = figure
        sub1 = subplot(1,2,1)
        pause(1)

        if(k == 1)
            %barDataToPlot = allAuc([3 4 5 6 1],[1 2])';% unfiltered   %[aucLinkItr ;aucLinkNoMuItr; aucUnLinkItr; aucUnLinkNoMuItr]';
            barDataToPlot = allAuc([3 4 5 6 7],[1 2])';% unfiltered   %[aucLinkItr ;aucLinkNoMuItr; aucUnLinkItr; aucUnLinkNoMuItr]';
            titleStr = ' All';
            figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_AUROC_allSites' '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
        else
            %barDataToPlot = allAuc([3 4 5 6 1],[3 4])';% Selc
            barDataToPlot = allAuc([3 4 5 6 7],[3 4])';% Selc
            titleStr = 'Selc';
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
        leg = legend('MPL', 'SL', 'MPLE-s_i', 'MPLE-s_{ij}', 'MPLE-R-s_i', 'location', 'NorthOutSide');
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
        leg = legend('MPL', 'SL', 'MPLE-s_i', 'MPLE-s_{ij}', 'MPLE-R-s_i', 'location', 'NorthOutSide');
        set(leg,'color','none');
        set(leg, 'Edgecolor','none');

        pause(1)
        % title(' ')


        if(saveFigs == 1)
            set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 8])%[0 0 19 6])% ,[0 0 8 6])
            set(gcf, 'renderer', 'painters');
            print([dirNameFigures figname], '-dpng','-r400')
            pause(0.5)
            %print(figname, '-depsc')
        end

    end





    %%
    % bar plot of ONLY MPLE 
    %----------------------------
    fig5 = figure
    pause(1)

    % barDataToPlot = allAuc(:,[1 2])';% unfiltered   %[aucLinkItr ;aucLinkNoMuItr; aucUnLinkItr; aucUnLinkNoMuItr]';
    % %barDataToPlot = allAuc(1,[3 4])';% Selc
    % barDataToPlot = [barDataToPlot barDataToPlot]';
    % bar(barDataToPlot)
    % 
    %myColorMap = color_scheme_npg([4 2 8],:);
    %color_scheme_npg
    %colormap([mapOranges(5,:); mapOranges(3,:); mapBlues(5,:); mapBlues(3,:)])
    %colormap([mapOranges(5,:); mapBlues(5,:); mapBlues(3,:)])
    % 


    for k = 1:2
        fig1 = figure
        sub1 = subplot(1,2,1)
        pause(1)

        if(k == 1)
            barDataToPlot = allAuc([5 6 7],[1 2])';% unfiltered   %[aucLinkItr ;aucLinkNoMuItr; aucUnLinkItr; aucUnLinkNoMuItr]';
            titleStr = ' All';
            figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_AUROC_allSites_OnlyMPLE_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
        else
            barDataToPlot = allAuc([5 6 7],[3 4])';% Selc
            titleStr = 'Selc';
            figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_AUROC_selcSites_OnlyMPLE_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
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

        %LabelBen = { 'Estimator (beneficial)'};
        %LabelDele = {'Estimator (deleterious)'};
        LabelBen = { 'Ben'};
        LabelDele = {'Del'};
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

        leg = legend('MPLE-s_i', 'MPLE-s_{ij}', 'MPLE-R-s_i', 'location', 'NorthOutSide');
        %leg = legend('MPLE-s_i', 'MPLE-s_{ij}', 'MPLE-All', 'location', 'NorthOutSide');

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
        leg = legend('MPLE-s_i', 'MPLE-s_{ij}', 'MPLE-R-s_i', 'location', 'NorthOutSide');
        set(leg,'color','none');
        set(leg, 'Edgecolor','none');

        pause(1)
        % title(' ')


        if(saveFigs == 1)
            set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 6 6])%[0 0 19 6])% ,[0 0 8 6])
            set(gcf, 'renderer', 'painters');
            print([dirNameFigures figname], '-dpng','-r400')
            pause(0.5)
            %print(figname, '-depsc')
        end

    end

    %%
    % make chord plots

    % myColorMap = color_scheme_npg([4 2 8],:);
    % 
    % alpha1 = 0.6;
    % 
    % myColorMap(1,:) = myColorMap(1,:)*(alpha1) + [1 1 1]*(1 - alpha1);
    % myColorMap(2,:) = myColorMap(2,:)*(alpha1) + [1 1 1]*(1 - alpha1);
    % myColorMap(3,:) = myColorMap(3,:)*(alpha1) + [1 1 1]*(1 - alpha1);

    alpha1 = 0.7;
    myColorMap2(1,:) = color_scheme_npg(4,:)*(alpha1) + [1 1 1]*(1 - alpha1);
    myColorMap2(2,:) = color_scheme_npg(2,:)*(alpha1) + [1 1 1]*(1 - alpha1);
    myColorMap2(3,:) = color_scheme_npg(8,:)*(alpha1) + [1 1 1]*(1 - alpha1);
    colormap(myColorMap2);

    % Create custom node labels
    myLabel = cell(5);
    for i = 1:5
      myLabel{i} = ['s_' num2str(i)];%['L_' num2str(i)];
    end

    perSiteSelctionEpiRd = round(perSiteSelctionEpi*10000)/10000;

    myColorMap2(4,:) = [0 0 0];
    figure
    circularGraph_noButtons(perSiteSelctionEpiRd*100,'Colormap',myColorMap2,'Label',myLabel);

    if(saveFigs == 1)
        %figname = ['/local/staff/ee/mssohail/Matlab Codes/H3N2 evolution/Fiftysites_3classPotts_Set' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_beneDele_all4methods_barChart'];
        figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_chordPlot_GT_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 4 4])% ,[0 0 8 6])
        %set(gcf, 'renderer', 'painters');
        print([dirNameFigures figname], '-dpng','-r400')
        %print(figname, '-depsc')
    end

    %%
    %selCoeffBoxPlot(allEstEpiItr, allEstItr, perSiteSelctionEpi, Lin)

    allEstEpiForBoxPlot = [];
    allClassEpiIndicatorForBoxPlot = [];
    allEstForBoxPlot = [];
    allClassIndicatorForBoxPlot = [];
    for l2 = 1:15
        thisSelcSiteEpiLogicInd = logical(allSelcSitesEpi(:,l2));
        allEstEpiForBoxPlot = [allEstEpiForBoxPlot ;allEstEpiItr(thisSelcSiteEpiLogicInd,l2)];
        allClassEpiIndicatorForBoxPlot = [allClassEpiIndicatorForBoxPlot ; repmat(l2,  sum(thisSelcSiteEpiLogicInd), 1)]; 
        if(sum(thisSelcSiteEpiLogicInd) == 0)
            allEstEpiForBoxPlot = [allEstEpiForBoxPlot ; -1];
            allClassEpiIndicatorForBoxPlot = [allClassEpiIndicatorForBoxPlot ; l2];
        end
        if(l2 < 6)
            thisSelcSiteLogicInd = logical(allSelcSites(:,l2));
            allEstForBoxPlot = [allEstForBoxPlot ;allEstItr(thisSelcSiteLogicInd,l2)];
            allClassIndicatorForBoxPlot = [allClassIndicatorForBoxPlot ; repmat(l2,  sum(thisSelcSiteLogicInd), 1)]; 
            if(sum(thisSelcSiteEpiLogicInd) == 0)
                allEstEpiForBoxPlot = [allEstEpiForBoxPlot ; -1];
                allClassEpiIndicatorForBoxPlot = [allClassEpiIndicatorForBoxPlot ; l2];
            end
        end
    end
    color_scheme = [repmat(color_scheme_npg(3,:), 5, 1); repmat(color_scheme_npg(5,:), 10, 1)];
    only1plot = 0;
    selCoeffBoxPlot_new(allEstEpiForBoxPlot, allClassEpiIndicatorForBoxPlot, allEstForBoxPlot, allClassIndicatorForBoxPlot, perSiteSelctionEpi, Lin, color_scheme, only1plot)

    if(saveFigs == 1)
        figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_boxPlot_Both' '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 15 12])%,[0 0 8 6.45])% ,[0 0 8 6])
        %set(gcf, 'renderer', 'painters');
        print([dirNameFigures figname], '-dpng','-r400')
    end

    only1plot = 1;
    selCoeffBoxPlot_new(allEstEpiForBoxPlot, allClassEpiIndicatorForBoxPlot, allEstForBoxPlot, allClassIndicatorForBoxPlot, perSiteSelctionEpi, Lin, color_scheme, only1plot)

    if(saveFigs == 1)
        figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_boxPlot_onlyMPLE' '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 4])%,[0 0 8 6.45])% ,[0 0 8 6])
        %set(gcf, 'renderer', 'painters');
        print([dirNameFigures figname], '-dpng','-r400')
    end


    %%% Warning: currently plots all site estimates, should only chose
    %%% selcSites for calculation
    % if(saveFigs == 1)
    %     figname = ['FiveSites_3class_Set' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_boxPlot_estEpiAve'];
    %     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 6 6])%,[0 0 8 6.45])% ,[0 0 8 6])
    %     set(gcf, 'renderer', 'painters');
    %     print(figname, '-dpng','-r400')
    % end

    %selCoeffBoxPlot(allEstEpiItr, allEstItr, perSiteSelctionEpi, Lin)

    if(saveFigs == 1)
    % figname = [mainDir chosenSlash 'Fig_5site_Set' num2str(thisSet) '_InitStrains' num2str(numStrainsInInitialPop) '_Selc'];
    % set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 13])%,[0 0 8 6.45])% ,[0 0 8 6])
    % set(gcf, 'renderer', 'painters');
    % print(figname, '-dpng','-r400')
    end

    %%%% does not work righ tnow, need to edit to only show estimate of Selc
    % values
    %selCoeffBoxPlot_Selc(allEstEpi, allEst, perSiteSelctionEpi, Lin, allSelcSitesEpi, allSelcSites)

    estSiteSelctionEpi = zeros(Lin, Lin);
    for l2 = 1:15
        thisSelcSiteLogicInd = logical(ones(numItr,1));%logical(allSelcSitesEpi(:,l2));
        if(l2 <= 5)
            estSiteSelctionEpi(l2,l2) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        elseif(l2 > 5 && l2 <= 9)
            estSiteSelctionEpi(1,l2-4) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        elseif(l2 > 9 && l2 <= 12)
            estSiteSelctionEpi(2,l2-7) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        elseif(l2 > 12 && l2 <= 14)
            estSiteSelctionEpi(3,l2-9) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        elseif(l2 > 14)
            estSiteSelctionEpi(4,l2-10) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
        end    
        allEstEpiAve(l2) = mean(allEstEpiItr(thisSelcSiteLogicInd,l2));
    end

    estSiteSelctionEpiRd = round(estSiteSelctionEpi*10000)/10000;
    estSiteSelctionEpiRd(isnan(estSiteSelctionEpiRd)) = 0;
    figure
    circularGraph_noButtons(estSiteSelctionEpiRd*100,'Colormap',myColorMap2,'Label',myLabel);

    if(saveFigs == 1)
        figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_chordPlot_estEpiAve' '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 4 4])%,[0 0 8 6.45])% ,[0 0 8 6])
        %set(gcf, 'renderer', 'painters');
        print([dirNameFigures figname], '-dpng','-r400')
    end



    %% make the chord diagram for a single itr

    chosenItrAll = [5:7]%13;

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
        circularGraph_noButtons(estSiteSelctionEpiChosenItrRd*100,'Colormap',myColorMap2,'Label',myLabel);

        if(saveFigs == 1)
            figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_chordPlot_estEpi_itr' num2str(chosenItr) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
            set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 3 3])%,[0 0 8 6.45])% ,[0 0 8 6])
            %set(gcf, 'renderer', 'painters');
            print([dirNameFigures figname], '-dpng','-r400')
        end
    end


    if(exportData2Excel)
        SelCoeff = selcCoeffNameAllItr;
        Estimates = estimatesAllItr;
        Method = methodsAllItr;
        Valid = validAllItr;
        TableOfSelecEst = table(SelCoeff, Estimates, Method, Valid);
        fileName_Table = [dirNameFigures chosenSlash 'TableEstimates_5Sites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.xlsx'];
        writetable(TableOfSelecEst, fileName_Table);
    end


    %%

    % two window bar plot

    % leftMargin = 0.2;
    % rightMargin = 0.02;
    % bottomMargin = 0.25;
    % topMargin = 0.12;
    % hgap = 0.03%0.09;
    % height2 = (1 - bottomMargin - topMargin);
    % height1 = height2;
    % width1 = (1 - leftMargin - rightMargin - hgap)/2;
    % width2 = width1;
    % 
    % 
    % fig1 = figure
    % ha(1) = axes('Units','normalized', ...
    %                 'Position',[leftMargin bottomMargin width1 height1], ...
    %                 'XTickLabel','', ...
    %                 'YTickLabel','');
    %             
    %             %
    %             
    %             
    % ha(2) = axes('Units','normalized', ...
    %                 'Position',[leftMargin+hgap+width1 bottomMargin width2 height2], ...
    %                 'XTickLabel','', ...
    %                 'YTickLabel','');
    % 
    %             
    % %LabelNumSiteSestSi = {'0.2', '0.4', '0.6', '0.8', '1'};
    % %LabelNumSiteSestSij = {'0.2', '0.4', '0.6', '0.8', '1'};
    % LabelNumSiteSestSi = {'0.2', '0.6', '1'};
    % LabelNumSiteSestSij = {'0.2', '0.6', '1'};
    % 
    % %h1 = histogram(numSelcSitesEpiItr, [0.5:1:16.5])
    % %h1 = histogram(numSelcSitesEpiItr_si, [0.5:1:16.5])
    % axes(ha(1))
    % h1 = histogram(numSelcSitesEpiItr_si/5, [0:.2:1.2])
    % h1.FaceColor = color_scheme_npg(3,:);
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
    %   'XTick', [0.3:0.4:1.2], ...
    %   'XTickLabel', LabelNumSiteSestSi, ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
    %   'LineWidth', 0.5)
    % ylabel('Count')
    % %xlabel('Number of accessible terms')
    % xlabel('Fraction of accessible terms')
    % xlabh = get(gca,'xlabel'); 
    % set(xlabh,'Units','normalized');
    % set(xlabh,'position',get(xlabh,'position') + [0.5 -0.1 0]);
    % ylim([0 1000])
    % axis([0 1.2 0 numItr])
    % 
    % 
    % axes(ha(2))
    % h2 = histogram(numSelcSitesEpiItr_sij/10, [0:.1:1.1])
    % h2.FaceColor = color_scheme_npg(5,:);
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
    % 
    % 
    % if(saveFigs == 1)
    %     figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_histNumTermsEst_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused)];
    %     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 3.5])%[0 0 6 4])%,[0 0 8 6.45])% ,[0 0 8 6])
    %     %set(gcf, 'renderer', 'painters');
    %     print([dirNameFigures figname], '-dpng','-r400')
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
    h1 = histogram(numSelcSitesEpiItr_si/5, [0:.2:1.2])
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
    axis([0 1.2 0 numItr])

    hold on
    %axes(ha(2))
    h2 = histogram(numSelcSitesEpiItr_sij/10, [0:.1:1.1])
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
        figname = ['FiveSites_Set' num2str(thisSet) '_numStrains' num2str(numStrainsInInitialPop) '_histNumTermsEst_ng' num2str(ng) '_dT' num2str(dT) '_Tused' num2str(Tused) '_reg1' regStr1  '_reg2' regStr2 '.png'];
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 5 3.5])%[0 0 6 4])%,[0 0 8 6.45])% ,[0 0 8 6])
        %set(gcf, 'renderer', 'painters');
        print([dirNameFigures figname], '-dpng','-r400')
    end



    %%% Saving data to send to John
    % dataToWrite = [sigmaEstOutLink_AllItr' trueClassROC_AllItr'];
    % 
    % dirNameGenData = [dirNameScriptFile '/FigData/'];
    % textFileNameToSave = 'Data_Fig2.csv';
    % csvwrite([dirNameGenData textFileNameToSave], dataToWrite);


    % %% Saving data to send to John
    %  dataToWrite = [sigmaEstOutLink_AllItr' trueClassROC_AllItr'];
    %  
    %  dirNameGenData = [dirNameScriptFile '/FigData/'];
    %  textFileNameToSave = 'Data_Fig_S14_c.csv';
    %  csvwrite([dirNameGenData textFileNameToSave], dataToWrite);

    %   dataToWrite = [sigmaEstOutUnLinkNoMu_AllItr' trueClassROC_AllItr'];
    %  
    %  dirNameGenData = [dirNameScriptFile '/FigData/'];
    %  textFileNameToSave = 'Data_Fig_S9_b.csv';
    %  csvwrite([dirNameGenData textFileNameToSave], dataToWrite);

    %%
    temp51 = strfind(fileName, 'itr');

    fileNamePlotDataSave = ['PlotData_' fileName(1:temp51+2) '1_' fileName(temp51+3:end)];
    save([dirNameFigures fileNamePlotDataSave])


end


%%

diffMtx = abs(sigmaEstOutLinkItr - sigmaEstOutLinkEpiItr(:,1:5));
diffVec = diffMtx(:);

figure
range2 = 0:0.0005:0.1;
[x,y]=ksdensity((diffVec),range2);

mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
%color_scheme3 = brewermap(100,'Greys');
color_scheme3 = brewermap(100,'Greens');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha2) + 255*white_scheme*(1-alpha2))/255;

blue = color_scheme21(50,:);
area(y-y(find(x==max(x))),(x),'FaceColor',blue,'FaceAlpha',0.7,'EdgeColor',blue,'LineWidth',1)
    