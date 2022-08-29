%% plots 2 of the  same meterics as previous but more plots more sets in same figure


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

thisSetAll = [418951001 418851001 418751001 418651001 418551001 418451001 418351001 418251001];
numStrainsInInitialPop = 10; % 5, 10, 30
repAll = [1 5 10 20]
lastItr = 1000;% for rep

fileNameContainingDirPath = 'dirNames.txt';
postFiltering = 1; % 0 : pre filtering (out put of AnalysisMPL_filt)
                   % 1 : post filtering
similarityThresh = 0.1; % threshold that sets allowed difference between two numbers to be still declared as same value, ie., 2 = 2.01
noiseThresh = 0.0;

dirNameScriptFile = pwd;    
if(ispc)
    chosenSlash = '\';
elseif(isunix)
    chosenSlash = '/';
else
    display('Error: system si not unix and not PC...')
    pause
end

dataPointCount = 1;
textCount = 1;
plotCounter = 1;
for i = 1:length(repAll)
    
    thisRep = repAll(i);
    for j = 1:length(thisSetAll)
        
        thisSet = thisSetAll(j); % 5, 10, 30
        
        getSysParam;
        Tstart = 1;%31;
        Tused = 100;%300;
        ng = 100;
        dT = 10;

        numParam = Lin*(Lin+1)/2;
        totalEpiTerms = Lin*(Lin-1)/2;
        step = 0.002;
        edges = [-10 -0.3:step:0.3 10] + step/2;


        %Tend = Tused + Tstart - 1;
        Tend = Tused + Tstart;

        actualT = T/1000;
        posSelc = selVal/2/Nin;
        negSelc = delSelVal/2/Nin;
        lineCol = {'r.-', 'b.-', 'k.-', 'g.-'};

        dataFilesWith2Regs = 1;

        regStr1 = '1';
        regStr2 = '1';
        useStrong = 0;

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


        if(thisRep == 1)
            numItr = lastItr;%100;

            
            numTps = Tused/dT;


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
            sigmaEstOutLinkEpiWithR_AllItr = [];
            sigmaEstOutLinkEpiWithRSelc_AllItr = [];
            sigmaEstOutLink_AllItr = [];
            sigmaEstOutLinkSelc_AllItr = [];
            sigmaEstOutUnLink_AllItr = [];
            sigmaEstOutUnLinkSelc_AllItr = [];
            perSiteSelction_AllItr  = [];
            perSiteSelctionAllEpiTerms_AllItr  = [];

            siteToExculeCuzOfLimitedPolyTps_AllItr = [];
            aucLinkEpiItr = -1*ones(numItr,2);
            aucLinkNoMuEpiItr = -1*ones(numItr,2);
            aucLinkEpiWithRItr = -1*ones(numItr,2);
            aucLinkItr = -1*ones(numItr,2);
            aucUnLinkItr = -1*ones(numItr,2);
            aucLinkEpiSelcItr = -1*ones(numItr,2);
            aucLinkNoMuEpiSelcItr = -1*ones(numItr,2);
            aucLinkEpiWithRSelcItr = -1*ones(numItr,2);
            aucLinkSelcItr = -1*ones(numItr,2);
            aucUnLinkSelcItr = -1*ones(numItr,2);

            numSelcSitesItr = zeros(1, numItr);
            numSelcSitesEpiItr = zeros(1, numItr);
            numSelcSitesEpiItr_si = zeros(1, numItr);
            numSelcSitesEpiItr_sij = zeros(1, numItr);

%             allEstEpiItr = zeros(numItr, numParam);
%             allEstNoMuEpiItr = zeros(numItr, numParam);
%             allEstEpiWithRItr = zeros(numItr, numParam);
%             allEstItr = zeros(numItr, Lin);

            allSelcSites = zeros(numItr, Lin);
            allSelcSitesEpi = zeros(numItr, numParam);
            allaccessibleAsSumSitesEpi = zeros(numItr, numParam);

            numTerms_LinkEpi_cluster_itr  = -1*ones(numItr,1);

            xCount = 1;

            aucLinkEpi_only_si_Itr = -1*ones(numItr, 2);
            aucLinkEpi_only_sij_Itr = -1*ones(numItr, 2);

            aucLinkEpi_only_si_SelcItr = -1*ones(numItr, 2);
            aucLinkEpi_only_sij_SelcItr = -1*ones(numItr, 2);


            aucLinkEpiWithR_only_si_Itr = -1*ones(numItr, 2);
            aucLinkEpiWithR_only_sij_Itr = -1*ones(numItr, 2);

            aucLinkEpiWithR_only_si_SelcItr = -1*ones(numItr, 2);
            aucLinkEpiWithR_only_sij_SelcItr = -1*ones(numItr, 2);


            nonSynSites = [];
            synSites = [];
            nonSynSitesEpi = [];
            synSitesEpi = [];

            for itr = 1:numItr%[1:50 52:100]%1:numItr
                itr
                selcSites = [];
                if(recombination == 1)
                    if(Tstart == 1)
                        fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
                    else        
                        fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
                    end
                    if(useStrong == 1)
                        fileName = [fileName(1:end-4) '_useStrong.mat'];
                    end
                    load([dirNameAnalysis fileName], 'Lin', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
                         'sigmaEstOutLink', 'sigmaEstOutLinkEpi', 'sigmaEstOutLinkEpiWithR', 'sigmaEstOutLinkNoMuEpi',...
                         'sigmaEstOutUnLink', 'sigmaEstOutLinkEpiNoE', ...
                         'q', 'q11', 'priorConst', 'timeWholeCodeMPL',...'perSiteSelctionSelc', 'selcSites', ... 'sigmaEstOutLinkSelc', 'sigmaEstOutLinkNoMuSelc','sigmaEstOutUnLinkSelc', 'sigmaEstOutUnLinkNoMuSelc', ...
                         'deltaLogLikeli', 'deltaLogLikeli50', 'deltaLogLikeli100',...
                         'deltaLogLikeliMPL', 'deltaLogLikeliSL', 'Hmat');

                else
                    if(Tstart == 1)
                        fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
                    else        
                        fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(itr) '.mat'];
                    end
                    if(useStrong == 1)
                        fileName = [fileName(1:end-4) '_useStrong.mat'];
                    end
                    load([dirNameAnalysis fileName], 'Lin', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
                         'sigmaEstOutLink', 'sigmaEstOutLinkEpi', 'sigmaEstOutLinkNoMuEpi',...
                         'sigmaEstOutUnLink', 'sigmaEstOutLinkEpiNoE', ...
                         'q', 'q11', 'priorConst', 'timeWholeCodeMPL',...'perSiteSelctionSelc', 'selcSites', ... 'sigmaEstOutLinkSelc', 'sigmaEstOutLinkNoMuSelc','sigmaEstOutUnLinkSelc', 'sigmaEstOutUnLinkNoMuSelc', ...
                         'deltaLogLikeli', 'deltaLogLikeli50', 'deltaLogLikeli100',...
                         'deltaLogLikeliMPL', 'deltaLogLikeliSL', 'Hmat');
                     sigmaEstOutLinkEpiWithR = zeros(length(sigmaEstOutLinkEpi), 1);
                end

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
                sitesWithZeroFreq = find(sum(Hmat) == 0);
                accessibleAsSumSites = setdiff(illCondSites, sitesWithZeroFreq);

                nonSynSites = [nonSynSites; sigmaEstOutLink(perSiteSelction ~= 0)];
                synSites = [synSites; sigmaEstOutLink(perSiteSelction == 0)];

                sigmaEstOutLinkEpiOnlyEpi = sigmaEstOutLinkEpi(1:Lin);
                nonSynSitesEpi = [nonSynSitesEpi; sigmaEstOutLinkEpiOnlyEpi(perSiteSelction ~= 0)];
                synSitesEpi = [synSitesEpi; sigmaEstOutLinkEpiOnlyEpi(perSiteSelction == 0)];

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
                    sigmaEstOutLinkEpiWithRSelc = sigmaEstOutLinkEpiWithR(selcSitesEpi);
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
                [aucLinkEpiWithRItrTemp, aucLinkEpiWithRSelcItrTemp] = calcAUROC(perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc, sigmaEstOutLinkEpiWithR, sigmaEstOutLinkEpiWithRSelc, posOnlyItrTemp_AllEpiTerms, negOnlyItrTemp_AllEpiTerms, posOnlySelcItrTemp_AllEpiTerms, negOnlySelcItrTemp_AllEpiTerms);

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


                % --- calculations for filtered ---

                [validItrTemp, xCellEntriesThisItrTemp, ...
                    xCellEntryTemp, yCellEntryTemp, aucCellEntryTemp] = calcAUROC_coord(neutralLimit, ...
                                                                     perSiteSelctionAllEpiTerms, selcSitesEpi_only_si, sigmaEstOutLinkEpi);

                [validItrTempWithR, xCellEntriesThisItrTempWithR, ...
                    xCellEntryTempWithR, yCellEntryTempWithR, aucCellEntryTempWithR] = calcAUROC_coord(neutralLimit, ...
                                                                     perSiteSelctionAllEpiTerms, selcSitesEpi_only_si, sigmaEstOutLinkEpiWithR);

%                 validItr(itr) = validItrTemp;
%                 xCellEntriesThisItr(itr) = xCellEntriesThisItrTemp;
%                 xCell{itr} = xCellEntryTemp;
%                 yCell{itr} = yCellEntryTemp;
%                 aucCell{itr} = aucCellEntryTemp;
% 
% 
%                 validItrWithR(itr) = validItrTempWithR;
%                 xCellEntriesThisItrWithR(itr) = xCellEntriesThisItrTempWithR;
%                 xCellWithR{itr} = xCellEntryTempWithR;
%                 yCellWithR{itr} = yCellEntryTempWithR;
%                 aucCellWithR{itr} = aucCellEntryTempWithR;


                zeroThresh = 0.001;%0.00025;
                naeReg = 0.001;%0.00025;

                % Save all variables for this itr
                %--------------------------------------

                % Save AUROC unfiltered
                aucLinkEpiItr(itr,:) = aucLinkEpiItrTemp;
                aucLinkNoMuEpiItr(itr,:) = aucLinkNoMuEpiItrTemp;
                aucLinkEpiWithRItr(itr,:) = aucLinkEpiWithRItrTemp;
                aucLinkItr(itr,:) = aucLinkItrTemp;
                aucUnLinkItr(itr,:) = aucUnLinkItrTemp;

                % Save AUROC filtered
                aucLinkEpiSelcItr(itr,:) = aucLinkEpiSelcItrTemp;
                aucLinkNoMuEpiSelcItr(itr,:) = aucLinkNoMuEpiSelcItrTemp;
                aucLinkEpiWithRSelcItr(itr,:) = aucLinkEpiWithRSelcItrTemp;
                aucLinkSelcItr(itr,:) = aucLinkSelcItrTemp;
                aucUnLinkSelcItr(itr,:) = aucUnLinkSelcItrTemp;

                aucLinkEpi_only_si_Itr(itr,:) = aucLinkEpiItrTemp_only_si;
                aucLinkEpi_only_sij_Itr(itr,:) = aucLinkEpiItrTemp_only_sij;

                aucLinkEpi_only_si_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_si;
                aucLinkEpi_only_sij_SelcItr(itr,:) = aucLinkEpiSelcItrTemp_only_sij;


                aucLinkEpiWithR_only_si_Itr(itr,:) = aucLinkEpiWithRItrTemp_only_si;
                aucLinkEpiWithR_only_sij_Itr(itr,:) = aucLinkEpiWithRItrTemp_only_sij;

                aucLinkEpiWithR_only_si_SelcItr(itr,:) = aucLinkEpiWithRSelcItrTemp_only_si;
                aucLinkEpiWithR_only_sij_SelcItr(itr,:) = aucLinkEpiWithRSelcItrTemp_only_sij;

                % Save estimated SC
                allEstEpiItr(itr,:) = sigmaEstOutLinkEpi';
                allEstNoMuEpiItr(itr,:) = sigmaEstOutLinkNoMuEpi';
                allEstEpiWithRItr(itr,:) = sigmaEstOutLinkEpiWithR';
                allEstItr(itr,:) = sigmaEstOutLink';    

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
            mean(aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_Itr(aucLinkEpi_only_sij_Itr(:,2) ~= -1,2)) mean(aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpi_only_sij_SelcItr(aucLinkEpi_only_sij_SelcItr(:,2) ~= -1,2));
            mean(aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_Itr(aucLinkEpiWithR_only_si_Itr(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,2) ~= -1,2));
            mean(aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_Itr(aucLinkEpiWithR_only_sij_Itr(:,2) ~= -1,2)) mean(aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,1) ~= -1,1)) mean(aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,2) ~= -1,2))]
            disp('-------------------------------')
            
            allAuc_ToPlot = allAuc;
            
            numSelcSitesItr_ToPlot(i,j) = mean(numSelcSitesItr);
            numSelcSitesEpiItr_si_ToPlot(i,j) = mean(numSelcSitesEpiItr_si);
            numSelcSitesEpiItr_sij_ToPlot(i,j) = mean(numSelcSitesEpiItr_sij);
            
            SE_numSelcSitesEpiItr(i,j) = std(numSelcSitesEpiItr(numSelcSitesEpiItr ~= 0)/190)/sqrt(sum(numSelcSitesEpiItr ~= 0));
            SE_numSelcSitesEpiItr_si(i,j) = std(numSelcSitesEpiItr_si(numSelcSitesEpiItr_si ~= 0)/Lin)/sqrt(sum(numSelcSitesEpiItr_si ~= 0));
            SE_numSelcSitesEpiItr_sij(i,j) = std(numSelcSitesEpiItr_sij(numSelcSitesEpiItr_sij ~= 0)/totalEpiTerms)/sqrt(sum(numSelcSitesEpiItr_sij ~= 0));
            
            thisGroup = 1;
            SE_EpiWithR_only_si_SelcBen(i,j) = std((aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr(:,thisGroup) ~= -1));
            SE_EpiWithR_only_sij_SelcBen(i,j) = std((aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr(:,thisGroup) ~= -1));
            SE_only_si_SelcBen(i,j) = std((aucLinkItr(aucLinkItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr(:,thisGroup) ~= -1));
            thisGroup = 2;
            SE_EpiWithR_only_si_SelcDel(i,j) = std((aucLinkEpiWithR_only_si_SelcItr(aucLinkEpiWithR_only_si_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr(:,thisGroup) ~= -1));
            SE_EpiWithR_only_sij_SelcDel(i,j) = std((aucLinkEpiWithR_only_sij_SelcItr(aucLinkEpiWithR_only_sij_SelcItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr(:,thisGroup) ~= -1));
            SE_only_si_SelcDel(i,j) = std((aucLinkItr(aucLinkItr(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr(:,thisGroup) ~= -1));

            
            
            mean_aucLinkEpiWithR_si_SelecBen(i,j) = allAuc(7,3);
            mean_aucLinkEpiWithR_si_SelecDel(i,j) = allAuc(7,4);
            mean_aucLinkEpiWithR_sij_SelecBen(i,j) = allAuc(8,3);
            mean_aucLinkEpiWithR_sij_SelecDel(i,j) = allAuc(8,4);
            
            mean_aucLink_si_SelecBen(i,j) = allAuc(3,3);
            mean_aucLink_si_SelecDel(i,j) = allAuc(3,4);
            
            dataPointCount = dataPointCount + 1;
            
        else
            
            thisItrStart = 1;
            thisItrEnd = thisRep;%5;%10;%4%10;%20;
            numRepComb = thisItrEnd - thisItrStart + 1;
            numItr = lastItr - thisItrStart + 1;

            numItr_rep = numItr/numRepComb;
      
            xCount = 1;
            
            numSelcSitesItr_rep = zeros(1, numItr/thisItrEnd);
            numSelcSitesEpiItr_rep = zeros(1, numItr/thisItrEnd);
            numSelcSitesEpiItr_si_rep = zeros(1, numItr/thisItrEnd);
            numSelcSitesEpiItr_sij_rep = zeros(1, numItr/thisItrEnd);

            
            aucLinkEpiItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkNoMuEpiItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkEpiWithRItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucUnLinkItr_rep = -1*ones(numItr/thisItrEnd, 2);

            aucLinkEpiSelcItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkNoMuEpiSelcItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkEpiWithRSelcItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkSelcItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucUnLinkSelcItr_rep = -1*ones(numItr/thisItrEnd, 2);

            aucLinkEpi_only_si_Itr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkEpi_only_sij_Itr_rep = -1*ones(numItr/thisItrEnd, 2);

            aucLinkEpi_only_si_SelcItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkEpi_only_sij_SelcItr_rep = -1*ones(numItr/thisItrEnd, 2);

            aucLinkEpiWithR_only_si_Itr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkEpiWithR_only_sij_Itr_rep = -1*ones(numItr/thisItrEnd, 2);

            aucLinkEpiWithR_only_si_SelcItr_rep = -1*ones(numItr/thisItrEnd, 2);
            aucLinkEpiWithR_only_sij_SelcItr_rep = -1*ones(numItr/thisItrEnd, 2);
                    
            nonSynSites_rep = [];
            synSites_rep = [];
            nonSynSitesEpi_rep = [];
            synSitesEpi_rep = [];

            %-------------------------------------------------------------------------------------------    
            % 1.2 Load data (combining replicates)
            
            for itr = 1:numItr
                if(rem(itr, thisItrEnd) == 0)    
                    itr

                    if(recombination == 1)
                        if(Tstart == 1)
                            fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff'  '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str((xCount-1)*numRepComb +thisItrStart) '_' num2str(xCount*numRepComb) '.mat'];
                        else        
                            fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str((xCount-1)*numRepComb +thisItrStart) '_' num2str(xCount*numRepComb) '.mat'];
                        end
                        if(useStrong == 1)
                            fileName = [fileName(1:end-4) '_useStrong.mat'];
                        end
                        load([dirNameAnalysis fileName], 'Lin', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
                             'sigmaEstOutLinkRep', 'sigmaEstOutLinkEpiRep', 'sigmaEstOutLinkEpiWithRRep', 'sigmaEstOutLinkNoMuEpiRep',...
                             'sigmaEstOutUnLinkRep', 'sigmaEstOutLinkEpiNoERep', 'HmatDTAll');
                    else
                        if(Tstart == 1)
                            fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff'  '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str((xCount-1)*numRepComb +thisItrStart) '_' num2str(xCount*numRepComb) '.mat'];
                        else        
                            fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' num2str(actualT) 'k_' num2str(itr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff' '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str((xCount-1)*numRepComb +thisItrStart) '_' num2str(xCount*numRepComb) '.mat'];
                        end
                        if(useStrong == 1)
                            fileName = [fileName(1:end-4) '_useStrong.mat'];
                        end
                        load([dirNameAnalysis fileName], 'Lin', 'perSiteSelction', 'perSiteSelctionEpi', 'muVal', ...
                             'sigmaEstOutLinkRep', 'sigmaEstOutLinkEpiRep', 'sigmaEstOutLinkNoMuEpiRep',...
                             'sigmaEstOutUnLinkRep', 'sigmaEstOutLinkEpiNoERep', 'HmatDTAll');
                         sigmaEstOutLinkEpiWithRRep = zeros(length(sigmaEstOutLinkEpiRep), 1);
                    end

                    sigmaEstOutLink_rep = sigmaEstOutLinkRep;
                    sigmaEstOutLinkEpi_rep = sigmaEstOutLinkEpiRep;
                    sigmaEstOutLinkNoMuEpi_rep = sigmaEstOutLinkNoMuEpiRep;
                    sigmaEstOutUnLink_rep = sigmaEstOutUnLinkRep;
                    sigmaEstOutLinkEpiNoE_rep = sigmaEstOutLinkEpiNoERep;
                    sigmaEstOutLinkEpiWithR_rep = sigmaEstOutLinkEpiWithRRep;
                    
                    perSiteSelctionAllEpiTerms = [];
                    for l = 1:Lin
                        perSiteSelctionAllEpiTerms = [perSiteSelctionAllEpiTerms perSiteSelctionEpi(l,l+1:end)];
                    end

                    perSiteSelctionAllEpiTerms = [perSiteSelction perSiteSelctionAllEpiTerms];

                    %==================================================================

                    [selcSitesEpi_rep, wellCondSites_rep, illCondSites_rep, cluster_rep] = findIndColsOfHmat(HmatDTAll);
                    %[selcSitesEpi_rep_f, wellCondSites_rep_f, illCondSites_rep_f, cluster_rep_f] = findIndColsOfHmat_fast(HmatDTAll);
                    %selcSites_rep = sum(allSelcSites((xCount-1)*numRepComb + 1:xCount*numRepComb,:)) > 0;
                    selcSites_rep = findIndColsOfHmat(HmatDTAll(1:Lin,1:Lin));
                    sitesWithZeroFreq_rep = find(sum(HmatDTAll) == 0);
                    accessibleAsSumSites_rep = setdiff(illCondSites_rep, sitesWithZeroFreq_rep);

                    numSelcSitesItr_rep(xCount) = sum(selcSites_rep);
                    numSelcSitesEpiItr_rep(xCount) = sum(selcSitesEpi_rep);
                    numSelcSitesEpiItr_si_rep(xCount) = sum(selcSitesEpi_rep(1:Lin));
                    numSelcSitesEpiItr_sij_rep(xCount) = sum(selcSitesEpi_rep(Lin+1:end));




                    perSiteSelctionSelc_rep = perSiteSelction(selcSites_rep);
                    perSiteSelctionAllEpiTermsSelc_rep = perSiteSelctionAllEpiTerms(selcSitesEpi_rep);

                    sigmaEstOutLinkEpiSelc_rep = sigmaEstOutLinkEpi_rep(selcSitesEpi_rep);
                    sigmaEstOutLinkEpiWithRSelc_rep = sigmaEstOutLinkEpiWithR_rep(selcSitesEpi_rep);
                    sigmaEstOutLinkNoMuEpiSelc_rep = sigmaEstOutLinkNoMuEpi_rep(selcSitesEpi_rep);
                    sigmaEstOutLinkEpiNoESelc_rep = sigmaEstOutLinkEpiNoE_rep(selcSitesEpi_rep);
                    sigmaEstOutLinkSelc_rep = sigmaEstOutLink_rep(selcSites_rep);
                    sigmaEstOutUnLinkSelc_rep = sigmaEstOutUnLink_rep(selcSites_rep);


                    neutralLimit = 0.005;

                    [aucLinkItrTemp_rep, aucLinkSelcItrTemp_rep, aucUnLinkItrTemp_rep, aucUnLinkSelcItrTemp_rep, ...
                    aucLinkEpiItrTemp_rep, aucLinkEpiSelcItrTemp_rep, aucLinkNoMuEpiItrTemp_rep, aucLinkNoMuEpiSelcItrTemp_rep, ...    
                    aucLinkEpiNoEItrTemp_rep, aucLinkEpiNoESelcItrTemp_rep, aucLinkEpiItrTemp_only_si_rep, aucLinkEpiSelcItrTemp_only_si_rep, ...
                    aucLinkEpiItrTemp_only_sij_rep, aucLinkEpiSelcItrTemp_only_sij_rep] = calcEpiRef_AUROC_new(perSiteSelction, perSiteSelctionSelc_rep, ...
                                        perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc_rep, selcSitesEpi_rep, ...
                                        sigmaEstOutLink_rep, sigmaEstOutLinkSelc_rep, sigmaEstOutUnLink_rep, sigmaEstOutUnLinkSelc_rep, ...
                                        sigmaEstOutLinkEpi_rep, sigmaEstOutLinkEpiSelc_rep, sigmaEstOutLinkNoMuEpi_rep, sigmaEstOutLinkNoMuEpiSelc_rep, ...
                                        sigmaEstOutLinkEpiNoE_rep, sigmaEstOutLinkEpiNoESelc_rep, neutralLimit);

                    [~, ~, ~, ~, ...
                    aucLinkEpiWithRItrTemp_rep, aucLinkEpiWithRSelcItrTemp_rep, ~, ~, ...    
                    ~, ~, aucLinkEpiWithRItrTemp_only_si_rep, aucLinkEpiWithRSelcItrTemp_only_si_rep, ...
                    aucLinkEpiWithRItrTemp_only_sij_rep, aucLinkEpiWithRSelcItrTemp_only_sij_rep] = calcEpiRef_AUROC_new(perSiteSelction, perSiteSelctionSelc_rep, ...
                                        perSiteSelctionAllEpiTerms, perSiteSelctionAllEpiTermsSelc_rep, selcSitesEpi_rep, ...
                                        sigmaEstOutLink_rep, sigmaEstOutLinkSelc_rep, sigmaEstOutUnLink_rep, sigmaEstOutUnLinkSelc_rep, ...
                                        sigmaEstOutLinkEpiWithR_rep, sigmaEstOutLinkEpiWithRSelc_rep, sigmaEstOutLinkNoMuEpi_rep, sigmaEstOutLinkNoMuEpiSelc_rep, ...
                                        sigmaEstOutLinkEpiNoE_rep, sigmaEstOutLinkEpiNoESelc_rep, neutralLimit);


                    zeroThresh = 0.001;%0.00025;
                    naeReg = 0.001;%0.00025;
            

%                     perSiteSelctionAllEpiTerms_cluster_rep = perSiteSelctionAllEpiTerms(wellCondSites_rep);
%                     sigmaEstOutLinkEpi_cluster_rep = sigmaEstOutLinkEpi_rep(wellCondSites_rep)';
%                     perSiteSelctionAllEpiTerms_clusterComb_rep = perSiteSelctionAllEpiTerms_cluster_rep;
%                     sigmaEstOutLinkEpi_clusterComb_rep = sigmaEstOutLinkEpi_cluster_rep;
%                     indTemp201_rep = [];
%                     num_clusterWellOnly_rep = length(sigmaEstOutLinkEpi_cluster_rep);
%                     for wt = 1:length(cluster_rep)
%                         indTemp101_rep = cluster_rep{wt};
%                         indTemp201_rep = [indTemp201_rep indTemp101_rep];
%                         perSiteSelctionAllEpiTerms_cluster_rep = [perSiteSelctionAllEpiTerms_cluster_rep sum(perSiteSelctionAllEpiTerms(indTemp101_rep))];
%                         sigmaEstOutLinkEpi_cluster_rep = [sigmaEstOutLinkEpi_cluster_rep sum(sigmaEstOutLinkEpi_rep(indTemp101_rep))];
%                     end
%                     clustInd_rep = zeros(1, length(sigmaEstOutLinkEpi_cluster_rep));
%                     temp65_rep = 1:num_clusterWellOnly_rep;
%                     clustInd_rep(temp65_rep) = 1;
%                     clustInd_rep = logical(clustInd_rep);



%                     perSiteSelctionAllEpiTerms_cluster_rep = round(100000*perSiteSelctionAllEpiTerms_cluster_rep)/100000;
%                     if(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep)) < zeroThresh)
%                         absErr_LinkEpi_cluster_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep - sigmaEstOutLinkEpi_cluster_rep))/sum(abs(perSiteSelctionAllEpiTerms_cluster_rep) + naeReg);%/length(perSiteSelctionAllEpiTerms_cluster_rep);
%                     else
%                         absErr_LinkEpi_cluster_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep - sigmaEstOutLinkEpi_cluster_rep))/sum(abs(perSiteSelctionAllEpiTerms_cluster_rep));%/length(perSiteSelctionAllEpiTerms_cluster_rep);
%                     end
%                     numTerms_LinkEpi_cluster_itr_rep(xCount) = length(perSiteSelctionAllEpiTerms_cluster_rep);
% 
%                     cluster_detail_cell_rep{itr,1} = perSiteSelctionAllEpiTerms_cluster_rep;
%                     cluster_detail_cell_rep{itr,2} = sigmaEstOutLinkEpi_cluster_rep;
% 
%                     if(sum(clustInd_rep) > 0)
%                         if(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep))) < zeroThresh)
%                             absErr_LinkEpi_clusterOnlyWell_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep) - sigmaEstOutLinkEpi_cluster_rep(clustInd_rep)))/(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep))) + naeReg);%/sum(clustInd_rep);
%                         else
%                             absErr_LinkEpi_clusterOnlyWell_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep) - sigmaEstOutLinkEpi_cluster_rep(clustInd_rep)))/sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(clustInd_rep)));%/sum(clustInd_rep);
%                         end
%                     end
% 
%                     if(sum(~clustInd_rep) > 0)
%                         if(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep))) < zeroThresh)
%                             absErr_LinkEpi_clusterOnlyAmb_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep) - sigmaEstOutLinkEpi_cluster_rep(~clustInd_rep)))/(sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep))) + naeReg);%/sum(~clustInd_rep);
%                         else
%                             absErr_LinkEpi_clusterOnlyAmb_itr_rep(xCount) = sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep) - sigmaEstOutLinkEpi_cluster_rep(~clustInd_rep)))/sum(abs(perSiteSelctionAllEpiTerms_cluster_rep(~clustInd_rep)));%/sum(~clustInd_rep);
%                         end
%                     end
% 
%                     if(isinf(absErr_LinkEpi_cluster_itr_rep(xCount)))
%                         pause
%                     end
% 
%                     if(absErr_LinkEpi_cluster_itr_rep(xCount) > 5)
%                         %pause
%                     end

                    % cluster     : well cond sites unique and ill cond sites divided into 
                    %               amb clusters, each cluster as 1 selCoeff
                    % clusterComb : all amb cluster bunched in 1 selCoeff



  

%                     allSelcSites_rep(xCount,:) = selcSites_rep;
%                     allSelcSitesEpi_rep(xCount,:) = selcSitesEpi_rep;

                    aucLinkEpiItr_rep(xCount,:) = aucLinkEpiItrTemp_rep;
                    aucLinkNoMuEpiItr_rep(xCount,:) = aucLinkNoMuEpiItrTemp_rep;
                    aucLinkEpiWithRItr_rep(xCount,:) = aucLinkEpiWithRItrTemp_rep;
                    aucLinkItr_rep(xCount,:) = aucLinkItrTemp_rep;
                    aucUnLinkItr_rep(xCount,:) = aucUnLinkItrTemp_rep;

                    % calculations for filtered
                    aucLinkEpiSelcItr_rep(xCount,:) = aucLinkEpiSelcItrTemp_rep;
                    aucLinkNoMuEpiSelcItr_rep(xCount,:) = aucLinkNoMuEpiSelcItrTemp_rep;
                    aucLinkEpiWithRSelcItr_rep(xCount,:) = aucLinkEpiWithRSelcItrTemp_rep;
                    aucLinkSelcItr_rep(xCount,:) = aucLinkSelcItrTemp_rep;
                    aucUnLinkSelcItr_rep(xCount,:) = aucUnLinkSelcItrTemp_rep;

                    aucLinkEpi_only_si_Itr_rep(xCount,:) = aucLinkEpiItrTemp_only_si_rep;
                    aucLinkEpi_only_sij_Itr_rep(xCount,:) = aucLinkEpiItrTemp_only_sij_rep;

                    aucLinkEpi_only_si_SelcItr_rep(xCount,:) = aucLinkEpiSelcItrTemp_only_si_rep;
                    aucLinkEpi_only_sij_SelcItr_rep(xCount,:) = aucLinkEpiSelcItrTemp_only_sij_rep;

                    aucLinkEpiWithR_only_si_Itr_rep(xCount,:) = aucLinkEpiWithRItrTemp_only_si_rep;
                    aucLinkEpiWithR_only_sij_Itr_rep(xCount,:) = aucLinkEpiWithRItrTemp_only_sij_rep;

                    aucLinkEpiWithR_only_si_SelcItr_rep(xCount,:) = aucLinkEpiWithRSelcItrTemp_only_si_rep;
                    aucLinkEpiWithR_only_sij_SelcItr_rep(xCount,:) = aucLinkEpiWithRSelcItrTemp_only_sij_rep;


                    condHmatWellItr_rep(itr) = cond(HmatDTAll(wellCondSites_rep,wellCondSites_rep));

                    nonSynSites_rep = [nonSynSites_rep; sigmaEstOutLinkRep(perSiteSelction ~= 0)];
                    synSites_rep = [synSites_rep; sigmaEstOutLinkRep(perSiteSelction == 0)];

                    sigmaEstOutLinkEpiOnlyEpiRep = sigmaEstOutLinkEpiRep(1:Lin);
                    nonSynSitesEpi_rep = [nonSynSitesEpi_rep; sigmaEstOutLinkEpiOnlyEpiRep(perSiteSelction ~= 0)];
                    synSitesEpi_rep = [synSitesEpi_rep; sigmaEstOutLinkEpiOnlyEpiRep(perSiteSelction == 0)];

%                     allEstEpiItr_rep(xCount,:) = sigmaEstOutLinkEpi_rep';
%                     allEstNoMuEpiItr_rep(xCount,:) = sigmaEstOutLinkNoMuEpi_rep';
%                     allEstEpiWithRItr_rep(xCount,:) = sigmaEstOutLinkEpiWithR_rep';
%                     allEstItr_rep(xCount,:) = sigmaEstOutLink_rep';



                            % --- calculations for filtered ---

                    % only s_i term sof Epi    
                    sitesEpi_only_si = logical([ones(1, Lin), zeros(1, Lin*(Lin-1)/2)]);
                    selcSitesEpi_only_si_rep = selcSitesEpi_rep & sitesEpi_only_si;

                    [validItrTemp_2, xCellEntriesThisItrTemp_2, ...
                        xCellEntryTemp_2, yCellEntryTemp_2, aucCellEntryTemp_2] = calcAUROC_coord(neutralLimit, ...
                                                                     perSiteSelctionAllEpiTerms, selcSitesEpi_only_si_rep, sigmaEstOutLinkEpi_rep);
                    [validItrTempWithR_2, xCellEntriesThisItrTempWithR_2, ...
                        xCellEntryTempWithR_2, yCellEntryTempWithR_2, aucCellEntryTempWithR_2] = calcAUROC_coord(neutralLimit, ...
                                                                     perSiteSelctionAllEpiTerms, selcSitesEpi_only_si_rep, sigmaEstOutLinkEpiWithR_rep);

%                     validItr_rep(xCount) = validItrTemp_2;
%                     xCellEntriesThisItr_rep(xCount) = xCellEntriesThisItrTemp_2;
%                     xCell_rep{xCount} = xCellEntryTemp_2;
%                     yCell_rep{xCount} = yCellEntryTemp_2;
%                     aucCell_rep{xCount} = aucCellEntryTemp_2;
% 
%                     validItrWithR_rep(xCount) = validItrTempWithR_2;
%                     xCellEntriesThisItrWithR_rep(xCount) = xCellEntriesThisItrTempWithR_2;
%                     xCellWithR_rep{xCount} = xCellEntryTempWithR_2;
%                     yCellWithR_rep{xCount} = yCellEntryTempWithR_2;
%                     aucCellWithR_rep{xCount} = aucCellEntryTempWithR_2;

                    xCount = xCount + 1;        
                end
            end
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

            whos aucLinkEpiWithR_only_si_Itr_rep
            pause(1)

            allAuc_ToPlot = allAuc_rep;

            numSelcSitesItr_ToPlot(i,j) = mean(numSelcSitesItr_rep);
            numSelcSitesEpiItr_si_ToPlot(i,j) = mean(numSelcSitesEpiItr_si_rep);
            numSelcSitesEpiItr_sij_ToPlot(i,j) = mean(numSelcSitesEpiItr_sij_rep);
            
            SE_numSelcSitesEpiItr(i,j) = std(numSelcSitesEpiItr_rep(numSelcSitesEpiItr_rep ~= 0)/190)/sqrt(sum(numSelcSitesEpiItr_rep ~= 0));
            SE_numSelcSitesEpiItr_si(i,j) = std(numSelcSitesEpiItr_si_rep(numSelcSitesEpiItr_si_rep ~= 0)/Lin)/sqrt(sum(numSelcSitesEpiItr_si_rep ~= 0));
            SE_numSelcSitesEpiItr_sij(i,j) = std(numSelcSitesEpiItr_sij_rep(numSelcSitesEpiItr_sij_rep ~= 0)/totalEpiTerms)/sqrt(sum(numSelcSitesEpiItr_sij_rep ~= 0));
            
            thisGroup = 1;
            SE_EpiWithR_only_si_SelcBen(i,j) = std((aucLinkEpiWithR_only_si_SelcItr_rep(aucLinkEpiWithR_only_si_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr_rep(:,thisGroup) ~= -1));
            SE_EpiWithR_only_sij_SelcBen(i,j) = std((aucLinkEpiWithR_only_sij_SelcItr_rep(aucLinkEpiWithR_only_sij_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr_rep(:,thisGroup) ~= -1));
            SE_only_si_SelcBen(i,j) = std((aucLinkItr_rep(aucLinkItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr_rep(:,thisGroup) ~= -1));
            thisGroup = 2;
            SE_EpiWithR_only_si_SelcDel(i,j) = std((aucLinkEpiWithR_only_si_SelcItr_rep(aucLinkEpiWithR_only_si_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_si_SelcItr_rep(:,thisGroup) ~= -1));
            SE_EpiWithR_only_sij_SelcDel(i,j) = std((aucLinkEpiWithR_only_sij_SelcItr_rep(aucLinkEpiWithR_only_sij_SelcItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkEpiWithR_only_sij_SelcItr_rep(:,thisGroup) ~= -1));
            SE_only_si_SelcDel(i,j) = std((aucLinkItr_rep(aucLinkItr_rep(:,thisGroup) ~= -1,thisGroup)))/sqrt(sum(aucLinkItr_rep(:,thisGroup) ~= -1));

            mean_aucLinkEpiWithR_si_SelecBen(i,j) = allAuc_rep(7,3);
            mean_aucLinkEpiWithR_si_SelecDel(i,j) = allAuc_rep(7,4);
            mean_aucLinkEpiWithR_sij_SelecBen(i,j) = allAuc_rep(8,3);
            mean_aucLinkEpiWithR_sij_SelecDel(i,j) = allAuc_rep(8,4);
            
            mean_aucLink_si_SelecBen(i,j) = allAuc_rep(3,3);
            mean_aucLink_si_SelecDel(i,j) = allAuc_rep(3,4);
            
            
            dataPointCount = dataPointCount + 1;
        end
    end
end
%%
save(['D:\New\MPL Epi\Figures\Sites20\PosDelSel\Sparsity_rep_20_sites_8sets.mat'])
%%
xDim = 20;
yDim = 10.5; %7;
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

leftMargin = 1.25;
rightMargin = 0.25;
bottomMargin = 1;%0.2;
topMargin = 1;
hgap1 = 0.2;
hgap2 = 1.05 + hgap1;
vgap1 = 0.2;
vgap2 = 1.25;

heightFree = (yDim - bottomMargin - topMargin - vgap1);
widthFree = (xDim - leftMargin - rightMargin - 5*hgap1 - hgap2);

height = heightFree/2;
width1 = widthFree/6;



ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin+2*(width1+hgap1)+hgap2 bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin+3*(width1+hgap1)+hgap2 bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(5) = axes('Units','centimeters', ...
                'Position',[leftMargin+4*(width1+hgap1)+hgap2 bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(6) = axes('Units','centimeters', ...
                'Position',[leftMargin+5*(width1+hgap1)+hgap2 bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
            
            
ha(7) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(8) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(9) = axes('Units','centimeters', ...
                'Position',[leftMargin+2*(width1+hgap1)+hgap2 bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(10) = axes('Units','centimeters', ...
                'Position',[leftMargin+3*(width1+hgap1)+hgap2 bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(11) = axes('Units','centimeters', ...
                'Position',[leftMargin+4*(width1+hgap1)+hgap2 bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(12) = axes('Units','centimeters', ...
                'Position',[leftMargin+5*(width1+hgap1)+hgap2 bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');




color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(100,'Greens');
color_scheme4 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme41 = (255*color_scheme4*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme_cell{1} = color_scheme11;
color_scheme_cell{2} = color_scheme11;
color_scheme_cell{3} = color_scheme31;
color_scheme_cell{4} = color_scheme41;

color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;  
columnInd = [0 0 1 2 3 4 0 0 1 2 3 4];
for k = [3:6 9:12]
    axes(ha(k))
    i = columnInd(k);
    if(k <= 6)
        plot(mean_aucLinkEpiWithR_si_SelecBen(:,i), '.-', 'color', color_scheme11(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        hold on
        errorbar(mean_aucLinkEpiWithR_si_SelecBen(:,i), SE_EpiWithR_only_si_SelcBen(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
        plot(mean_aucLinkEpiWithR_sij_SelecBen(:,i), '--', 'color', color_scheme11(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        errorbar(mean_aucLinkEpiWithR_sij_SelecBen(:,i), SE_EpiWithR_only_sij_SelcBen(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    else
        plot(mean_aucLinkEpiWithR_si_SelecDel(:,i), '-', 'color', color_scheme21(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        hold on
        errorbar(mean_aucLinkEpiWithR_si_SelecDel(:,i), SE_EpiWithR_only_si_SelcDel(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
        plot(mean_aucLinkEpiWithR_sij_SelecDel(:,i), '--', 'color', color_scheme21(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        errorbar(mean_aucLinkEpiWithR_sij_SelecDel(:,i), SE_EpiWithR_only_sij_SelcDel(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    end
    
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
      'YTick', 0.5:0.1:1, ...
      'XTick', [1:1:6 8], ...
      'XTickLabel'  , '', ...
      'YTickLabel'  , '', ...
      'LineWidth', 0.5)
    if((k == 3 || k == 9))
        set(gca, ...
        'YTickLabel', 0.5:0.1:1)
    end
    if(k > 6)
        set(gca, ...
            'XTickLabel', [1 5 10 20])
    end
    if(k == 3)
        ylabel('AUROC (beneficial)')
    end
    if(k == 9)
        ylabel('AUROC (deleterious)')
    end
    
    axis([0.5 4.5 0.58 1])
    %ylabel('AUROC (beneficial)')

end
xlabel('Number of replicates')
xlabh = get(gca,'xlabel'); 
set(xlabh,'Units','centimeters');
set(xlabh,'position',get(xlabh,'position') + [-4.25 0 0]);

for k = [1 2 7 8]
    axes(ha(k))
    
    plot(numSelcSitesEpiItr_si_ToPlot(:,i)/Lin, '.-', 'color', color_scheme41(80,:), 'MarkerFaceColor', color_scheme41(80,:), 'lineWidth', 1)
    hold on
    errorbar(numSelcSitesEpiItr_si_ToPlot(:,i)/Lin, SE_numSelcSitesEpiItr_si(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    plot(numSelcSitesEpiItr_sij_ToPlot(:,i)/totalEpiTerms, '--', 'color', color_scheme41(80,:), 'MarkerFaceColor', color_scheme41(80,:), 'lineWidth', 1)
    errorbar(numSelcSitesEpiItr_sij_ToPlot(:,i)/totalEpiTerms, SE_numSelcSitesEpiItr_sij(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)

    
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
      'YTick', 0:0.2:1, ...
      'XTick', [1:1:6 8], ...
      'XTickLabel'  , '', ...
      'YTickLabel'  , '', ...
      'LineWidth', 0.5)
    if((k == 1 || k == 7))
        set(gca, ...
        'YTickLabel', 0:0.2:1)
    end
    if(k > 6)
        set(gca, ...
            'XTickLabel', [1 5 10 20])
    end
    axis([0.5 4.5 -0.02 1])
    

end
xlabel('Number of replicates')
xlabh2 = get(gca,'xlabel'); 
set(xlabh2,'Units','centimeters');
set(xlabh2,'position',get(xlabh2,'position') + [-1.35 0 0]);

ylabel('Fraction of accessible terms')
ylabh2 = get(gca,'ylabel'); 
set(ylabh2,'Units','centimeters');
set(ylabh2,'position',get(ylabh2,'position') + [-(width1+hgap1+0.5) 2.7 0]);

dimDummy = [1 1 1 1];


la = annotation('line',[0 0.1], [0 0.1]);
la.Units = 'centimeter';
la.X = [leftMargin+0.5 leftMargin+0.5+0.6];
la.Y = [yDim-0.75 yDim-0.75];
la.LineWidth = 1;
la.Color = color_scheme41(80,:);

% lb = annotation('line',[0 0.1], [0 0.1]);
% lb.Units = 'centimeter';
% lb.X = [leftMargin+0.5+3 leftMargin+0.5+0.6+3];
% lb.Y = [yDim-0.75 yDim-0.75];
% lb.LineWidth = 1;
% lb.Color = color_scheme41(80,:);
% lb.LineStyle = '--';
 
tLeg1 = annotation('textbox',dimDummy,'String','selection coefficients','FitBoxToText','on');
tLeg1.Units = 'centimeter';
tLeg1.Position = [leftMargin+0.5+0.6 yDim-0.7 0.7 0.5];
tLeg1.LineStyle = 'none';
tLeg1.FontName = 'Arial';
tLeg1.FontSize = 8;      
        
% tLeg2 = annotation('textbox',dimDummy,'String','selection coefficients','FitBoxToText','on');
% tLeg2.Units = 'centimeter';
% tLeg2.Position = [leftMargin+0.5+0.9+2.6 yDim-0.7 0.7 0.5];
% tLeg2.LineStyle = 'none';
% tLeg2.FontName = 'Arial';
% tLeg2.FontSize = 8;      
        








ta = annotation('textbox',dimDummy,'String','A','FitBoxToText','on')
ta.Units = 'centimeter';
ta.Position = [-0.14 yDim-0.5 0.7 0.5];
ta.LineStyle = 'none';
ta.FontWeight = 'bold';
ta.FontName = 'Arial';
ta.FontSize = 12;

tb = annotation('textbox',dimDummy,'String','B','FitBoxToText','on')
tb.Units = 'centimeter';
tb.Position = [leftMargin+2*width1+hgap1 yDim-0.5 0.7 0.5];
tb.LineStyle = 'none';
tb.FontWeight = 'bold';
tb.FontName = 'Arial';
tb.FontSize = 12;

asd
%%

%%
xDim = 12;
yDim = 14; %7;
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

leftMargin = 1.5;
rightMargin = 0.25;
bottomMargin = 1;%0.2;
topMargin = 1;
hgap1 = 0.2;
hgap2 = 1.05 + hgap1;
vgap1 = 0.2;
vgap2 = 1.25;

heightFree = (yDim - bottomMargin - topMargin - 2*vgap1);
widthFree = (xDim - leftMargin - rightMargin - 3*hgap1);

height = heightFree/3;
width1 = widthFree/4;



ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+2*(height+vgap1) width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin+2*(height+vgap1) width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin+2*(width1+hgap1) bottomMargin+2*(height+vgap1) width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin+3*(width1+hgap1) bottomMargin+2*(height+vgap1) width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(5) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(6) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(7) = axes('Units','centimeters', ...
                'Position',[leftMargin+2*(width1+hgap1) bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(8) = axes('Units','centimeters', ...
                'Position',[leftMargin+3*(width1+hgap1) bottomMargin+height+vgap1 width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(9) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(10) = axes('Units','centimeters', ...
                'Position',[leftMargin+width1+hgap1 bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(11) = axes('Units','centimeters', ...
                'Position',[leftMargin+2*(width1+hgap1) bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
ha(12) = axes('Units','centimeters', ...
                'Position',[leftMargin+3*(width1+hgap1) bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
            


color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(100,'Greens');
color_scheme4 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme41 = (255*color_scheme4*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme_cell{1} = color_scheme11;
color_scheme_cell{2} = color_scheme11;
color_scheme_cell{3} = color_scheme31;
color_scheme_cell{4} = color_scheme41;

color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;  
columnInd = [0 0 0 0 1 2 3 4 1 2 3 4];
for k = [5:12]
    axes(ha(k))
    i = columnInd(k);
    if(k <= 8)
        plot(mean_aucLinkEpiWithR_si_SelecBen(:,i), '.-', 'color', color_scheme11(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        hold on
        errorbar(mean_aucLinkEpiWithR_si_SelecBen(:,i), SE_EpiWithR_only_si_SelcBen(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
        plot(mean_aucLinkEpiWithR_sij_SelecBen(:,i), '--', 'color', color_scheme11(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        errorbar(mean_aucLinkEpiWithR_sij_SelecBen(:,i), SE_EpiWithR_only_sij_SelcBen(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    else
        plot(mean_aucLinkEpiWithR_si_SelecDel(:,i), '-', 'color', color_scheme21(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        hold on
        errorbar(mean_aucLinkEpiWithR_si_SelecDel(:,i), SE_EpiWithR_only_si_SelcDel(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
        plot(mean_aucLinkEpiWithR_sij_SelecDel(:,i), '--', 'color', color_scheme21(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        errorbar(mean_aucLinkEpiWithR_sij_SelecDel(:,i), SE_EpiWithR_only_sij_SelcDel(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    end
    
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
      'YTick', 0.5:0.1:1, ...
      'XTick', [1:1:6 8], ...
      'XTickLabel'  , '', ...
      'YTickLabel'  , '', ...
      'LineWidth', 0.5)
    if((k == 5 || k == 9))
        set(gca, ...
        'YTickLabel', 0.5:0.1:1)
    end
    if(k > 8)
        set(gca, ...
            'XTickLabel', [1 5 10 20])
    end
    if(k == 5)
        ylabel({'AUROC', '(beneficial)'})
    end
    if(k == 9)
        ylabel({'AUROC', '(deleterious)'})
    end
    
    axis([0.5 4.5 0.58 1])
    %ylabel('AUROC (beneficial)')

end
xlabel('Number of replicates')
xlabh = get(gca,'xlabel'); 
set(xlabh,'Units','centimeters');
set(xlabh,'position',get(xlabh,'position') + [-(1*width1+1*hgap1+1.25) 0 0]);
sparsityVec = 100 - [95 85 75 65];
for i = [1:4]
    axes(ha(i))
    
    plot(numSelcSitesEpiItr_si_ToPlot(:,i)/Lin, '.-', 'color', color_scheme41(80,:), 'MarkerFaceColor', color_scheme41(80,:), 'lineWidth', 1)
    hold on
    errorbar(numSelcSitesEpiItr_si_ToPlot(:,i)/Lin, SE_numSelcSitesEpiItr_si(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    plot(numSelcSitesEpiItr_sij_ToPlot(:,i)/totalEpiTerms, '--', 'color', color_scheme41(80,:), 'MarkerFaceColor', color_scheme41(80,:), 'lineWidth', 1)
    errorbar(numSelcSitesEpiItr_sij_ToPlot(:,i)/totalEpiTerms, SE_numSelcSitesEpiItr_sij(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)

    
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
      'YTick', 0:0.2:1, ...
      'XTick', [1:1:6 8], ...
      'XTickLabel'  , '', ...
      'YTickLabel'  , '', ...
      'LineWidth', 0.5)
    if(i == 1)
        set(gca, ...
        'YTickLabel', 0:0.2:1)
    end
    title(['Sparsity: ' num2str(sparsityVec(i)) '%'])
    axis([0.5 4.5 -0.06 1])
end

ylabel({'Fraction of', 'accessible terms'})
ylabh2 = get(gca,'ylabel'); 
set(ylabh2,'Units','centimeters');
set(ylabh2,'position',get(ylabh2,'position') + [-(3*width1+2*hgap1+0.75) 0 0]);

dimDummy = [1 1 1 1];


la = annotation('line',[0 0.1], [0 0.1]);
la.Units = 'centimeter';
la.X = [leftMargin+0.4 leftMargin+0.4+0.6]+0.2;
la.Y = [yDim-0.25 yDim-0.25];
la.LineWidth = 1;
la.Color = color_scheme41(60,:);

lb = annotation('line',[0 0.1], [0 0.1]);
lb.Units = 'centimeter';
lb.X = [leftMargin+0.4 leftMargin+0.4+0.6]+5.8;
lb.Y = [yDim-0.25 yDim-0.25];
lb.LineWidth = 1;
lb.Color = color_scheme41(60,:);
lb.LineStyle = '--';
 
tLeg1 = annotation('textbox',dimDummy,'String','Selection coefficients','FitBoxToText','on');
tLeg1.Units = 'centimeter';
tLeg1.Position = [leftMargin+0.4+0.6+0.2 yDim-0.45 5 0.5];
tLeg1.LineStyle = 'none';
tLeg1.FontName = 'Arial';
tLeg1.FontSize = 8;      
        
tLeg2 = annotation('textbox',dimDummy,'String','Epistasis terms','FitBoxToText','on');
tLeg2.Units = 'centimeter';
tLeg2.Position = [leftMargin+0.4+0.6+5.8 yDim-0.45 5 0.5];
tLeg2.LineStyle = 'none';
tLeg2.FontName = 'Arial';
tLeg2.FontSize = 8;      


figureDir = [pwd '\Figures\'];
if(saveFigs == 1)
    figname = [figureDir 'Figure_20Sites_numAccParam_sparsity'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    
     set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
end

%%
% % plot 8-set
%%
xDim = 20;
yDim = 12; %7;
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

leftMargin = 1.5;
rightMargin = 0.25;
bottomMargin = 1;%0.2;
topMargin = 1;
hgap1 = 0.2;
hgap2 = 1.05 + hgap1;
vgap1 = 0.2;
vgap2 = 1.25;

heightFree = (yDim - bottomMargin - topMargin - 2*vgap1);
widthFree = (xDim - leftMargin - rightMargin - 7*hgap1);

height = heightFree/3;
width1 = widthFree/8;


for ind1 = 1:8
    ha(ind1) = axes('Units','centimeters', ...
                    'Position',[leftMargin+(ind1-1)*(width1+hgap1) bottomMargin+2*(height+vgap1) width1 height], ...
                    'XTickLabel','', ...
                    'YTickLabel','');
end

for ind1 = 9:16
    ha(ind1) = axes('Units','centimeters', ...
                    'Position',[leftMargin+(ind1 - 9)*(width1+hgap1) bottomMargin+height+vgap1 width1 height], ...
                    'XTickLabel','', ...
                    'YTickLabel','');
end

for ind1 = 17:24
    ha(ind1) = axes('Units','centimeters', ...
                    'Position',[leftMargin+(ind1 - 17)*(width1+hgap1) bottomMargin width1 height], ...
                    'XTickLabel','', ...
                    'YTickLabel','');
end


color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(100,'Greens');
color_scheme4 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme41 = (255*color_scheme4*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme_cell{1} = color_scheme11;
color_scheme_cell{2} = color_scheme11;
color_scheme_cell{3} = color_scheme31;
color_scheme_cell{4} = color_scheme41;

color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;  
columnInd = [zeros(1, 8) 1:8 1:8];
for k = [9:24]
    axes(ha(k))
    i = columnInd(k);
    if(k <= 16)
        plot(mean_aucLinkEpiWithR_si_SelecBen(:,i), '.-', 'color', color_scheme11(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        hold on
        errorbar(mean_aucLinkEpiWithR_si_SelecBen(:,i), SE_EpiWithR_only_si_SelcBen(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
        plot(mean_aucLinkEpiWithR_sij_SelecBen(:,i), '--', 'color', color_scheme11(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        errorbar(mean_aucLinkEpiWithR_sij_SelecBen(:,i), SE_EpiWithR_only_sij_SelcBen(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    else
        plot(mean_aucLinkEpiWithR_si_SelecDel(:,i), '-', 'color', color_scheme21(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        hold on
        errorbar(mean_aucLinkEpiWithR_si_SelecDel(:,i), SE_EpiWithR_only_si_SelcDel(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
        plot(mean_aucLinkEpiWithR_sij_SelecDel(:,i), '--', 'color', color_scheme21(80,:), 'MarkerFaceColor', color_scheme11(80,:), 'lineWidth', 1)
        errorbar(mean_aucLinkEpiWithR_sij_SelecDel(:,i), SE_EpiWithR_only_sij_SelcDel(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    end
    
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
      'YTick', 0.5:0.1:1, ...
      'XTick', [1:1:6 8], ...
      'XTickLabel'  , '', ...
      'YTickLabel'  , '', ...
      'LineWidth', 0.5)
    if((k == 9 || k == 17))
        set(gca, ...
        'YTickLabel', 0.5:0.1:1)
    end
    if(k > 16)
        set(gca, ...
            'XTickLabel', [1 5 10 20])
    end
    if(k == 9)
        ylabel({'AUROC', '(beneficial)'})
    end
    if(k == 17)
        ylabel({'AUROC', '(deleterious)'})
    end
    
    axis([0.5 4.5 0.58 1])
    %ylabel('AUROC (beneficial)')

end
xlabel('Number of replicates')
xlabh = get(gca,'xlabel'); 
set(xlabh,'Units','centimeters');
set(xlabh,'position',get(xlabh,'position') + [-(3*width1+3*hgap1+1.05) 0 0]);
sparsityVec = 100 - [95:-10:25];
for i = [1:8]
    axes(ha(i))
    
    plot(numSelcSitesEpiItr_si_ToPlot(:,i)/Lin, '.-', 'color', color_scheme41(80,:), 'MarkerFaceColor', color_scheme41(80,:), 'lineWidth', 1)
    hold on
    errorbar(numSelcSitesEpiItr_si_ToPlot(:,i)/Lin, SE_numSelcSitesEpiItr_si(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    plot(numSelcSitesEpiItr_sij_ToPlot(:,i)/totalEpiTerms, '--', 'color', color_scheme41(80,:), 'MarkerFaceColor', color_scheme41(80,:), 'lineWidth', 1)
    errorbar(numSelcSitesEpiItr_sij_ToPlot(:,i)/totalEpiTerms, SE_numSelcSitesEpiItr_sij(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)

    
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
      'YTick', 0:0.2:1, ...
      'XTick', [1:1:6 8], ...
      'XTickLabel'  , '', ...
      'YTickLabel'  , '', ...
      'LineWidth', 0.5)
    if(i == 1)
        set(gca, ...
        'YTickLabel', 0:0.2:1)
    end
    title(['Sparsity: ' num2str(sparsityVec(i)) '%'])
    axis([0.5 4.5 -0.06 1])
end

ylabel({'Fraction of', 'accessible terms'})
ylabh2 = get(gca,'ylabel'); 
set(ylabh2,'Units','centimeters');
set(ylabh2,'position',get(ylabh2,'position') + [-(7*width1+6*hgap1+0.75) 0 0]);

dimDummy = [1 1 1 1];


la = annotation('line',[0 0.1], [0 0.1]);
la.Units = 'centimeter';
la.X = [leftMargin+0.4 leftMargin+0.4+0.6]+0.2;
la.Y = [yDim-0.25 yDim-0.25];
la.LineWidth = 1;
la.Color = color_scheme41(60,:);

lb = annotation('line',[0 0.1], [0 0.1]);
lb.Units = 'centimeter';
lb.X = [leftMargin+0.4 leftMargin+0.4+0.6]+5.8;
lb.Y = [yDim-0.25 yDim-0.25];
lb.LineWidth = 1;
lb.Color = color_scheme41(60,:);
lb.LineStyle = '--';
 
tLeg1 = annotation('textbox',dimDummy,'String','Selection coefficients','FitBoxToText','on');
tLeg1.Units = 'centimeter';
tLeg1.Position = [leftMargin+0.4+0.6+0.2 yDim-0.45 5 0.5];
tLeg1.LineStyle = 'none';
tLeg1.FontName = 'Arial';
tLeg1.FontSize = 8;      
        
tLeg2 = annotation('textbox',dimDummy,'String','Epistasis terms','FitBoxToText','on');
tLeg2.Units = 'centimeter';
tLeg2.Position = [leftMargin+0.4+0.6+5.8 yDim-0.45 5 0.5];
tLeg2.LineStyle = 'none';
tLeg2.FontName = 'Arial';
tLeg2.FontSize = 8;      


figureDir = [pwd '\Figures\'];
if(saveFigs == 1)
    figname = [figureDir 'FigureSX_20Sites_numAccParam_sparsity_8set'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    
     set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
end


%% plot only number of accessible figure

xDim = 20;
yDim = 5.6; %7;
fig1 = figure('Units','centimeters', ...
                'Position', [1 1 xDim yDim])

leftMargin = 1.5;
rightMargin = 0.25;
bottomMargin = 1;%0.2;
topMargin = 1;
hgap1 = 0.2;
hgap2 = 1.05 + hgap1;
vgap1 = 0.2;
vgap2 = 1.25;

heightFree = (yDim - bottomMargin - topMargin);
widthFree = (xDim - leftMargin - rightMargin - 7*hgap1);

height = heightFree;
width1 = widthFree/8;


for ind1 = 1:8
    ha(ind1) = axes('Units','centimeters', ...
                    'Position',[leftMargin+(ind1-1)*(width1+hgap1) bottomMargin width1 height], ...
                    'XTickLabel','', ...
                    'YTickLabel','');
end


color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(100,'Greens');
color_scheme4 = brewermap(100,'Greys');
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme41 = (255*color_scheme4*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme_cell{1} = color_scheme11;
color_scheme_cell{2} = color_scheme11;
color_scheme_cell{3} = color_scheme31;
color_scheme_cell{4} = color_scheme41;

color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;  

sparsityVec = 100 - [95:-10:25];
for i = [1:8]
    axes(ha(i))
    
    plot(numSelcSitesEpiItr_si_ToPlot(:,i)/Lin, '.-', 'color', color_scheme41(80,:), 'MarkerFaceColor', color_scheme41(80,:), 'lineWidth', 1)
    hold on
    errorbar(numSelcSitesEpiItr_si_ToPlot(:,i)/Lin, SE_numSelcSitesEpiItr_si(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    plot(numSelcSitesEpiItr_sij_ToPlot(:,i)/totalEpiTerms, '--', 'color', color_scheme41(80,:), 'MarkerFaceColor', color_scheme41(80,:), 'lineWidth', 1)
    errorbar(numSelcSitesEpiItr_sij_ToPlot(:,i)/totalEpiTerms, SE_numSelcSitesEpiItr_sij(:,i), 'k', 'LineStyle', 'none', 'CapSize', 3)
    %plot(0:0.5:5, ones(1, length(0:0.5:5)), 'k:')

    
    set(gca, ...
      'Box'         , 'on'     , ...
      'TickDir'     , 'in'     , ...
      'TickLength'  , [.01 .01] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off'      , ...
      'XGrid'       , 'off'      , ...
      'XColor'      , [.1 .1 .1], ...
      'YColor'      , [.1 .1 .1], ...'LineWidth'   , 1      , ...
      'YTick', 0:0.2:1, ...
      'XTick', [1:1:6 8], ...
     'XTickLabel', [1 5 10 20], ...
      'YTickLabel'  , '', ...
      'LineWidth', 0.5)
    if(i == 1)
        set(gca, ...
        'YTickLabel', 0:0.2:1)
    end
    title(['Sparsity: ' num2str(sparsityVec(i)) '%'])
    axis([0.5 4.5 -0.06 1])
end

xlabel('Number of replicates')
xlabh = get(gca,'xlabel'); 
set(xlabh,'Units','centimeters');
set(xlabh,'position',get(xlabh,'position') + [-(3*width1+3*hgap1+1.05) 0 0]);

ylabel({'Fraction of', 'accessible terms'})
ylabh2 = get(gca,'ylabel'); 
set(ylabh2,'Units','centimeters');
set(ylabh2,'position',get(ylabh2,'position') + [-(7*width1+6*hgap1+0.75) 0 0]);

dimDummy = [1 1 1 1];


la = annotation('line',[0 0.1], [0 0.1]);
la.Units = 'centimeter';
la.X = [leftMargin+0.4 leftMargin+0.4+0.6]+0.2;
la.Y = [yDim-0.25 yDim-0.25];
la.LineWidth = 1;
la.Color = color_scheme41(60,:);

lb = annotation('line',[0 0.1], [0 0.1]);
lb.Units = 'centimeter';
lb.X = [leftMargin+0.4 leftMargin+0.4+0.6]+5.8;
lb.Y = [yDim-0.25 yDim-0.25];
lb.LineWidth = 1;
lb.Color = color_scheme41(60,:);
lb.LineStyle = '--';
 
tLeg1 = annotation('textbox',dimDummy,'String','Selection coefficients','FitBoxToText','on');
tLeg1.Units = 'centimeter';
tLeg1.Position = [leftMargin+0.4+0.6+0.2 yDim-0.45 5 0.5];
tLeg1.LineStyle = 'none';
tLeg1.FontName = 'Arial';
tLeg1.FontSize = 8;      
        
tLeg2 = annotation('textbox',dimDummy,'String','Epistasis terms','FitBoxToText','on');
tLeg2.Units = 'centimeter';
tLeg2.Position = [leftMargin+0.4+0.6+5.8 yDim-0.45 5 0.5];
tLeg2.LineStyle = 'none';
tLeg2.FontName = 'Arial';
tLeg2.FontSize = 8;      


figureDir = [pwd '\Figures\'];
if(saveFigs == 1)
    figname = [figureDir 'FigureS10_20Sites_numAccParam_sparsity_8set'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
    
     set(gcf, 'PaperSize', [xDim yDim])
    print(figname, '-dpdf', '-fillpage')
end

