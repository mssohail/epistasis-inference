%% Run Maximum Path Likelihood (MPLin) method on GT data
% v2 is faster version. It calculates the int cov matrix in step 2 and does
% not save cov matrix at each time point

% Use different reg for epi terms
%%

% this code :
%      1. loads data (in .mat or .dat format)
%      2. calculates one and two point frequencies
%      3. Find estimates using MPL

% Last updated 03-Jan-2021
% code for 5-site simulation Figure 4

% previous update 10-Jan 2018
% Author: M Saqib Sohail


% 
clc

clear all
close all

saveFile = 1;%1;
useStrongReg = 0; % if 1, uses 10x stronger reg for epi sites that are not accessible
numStrainsInInitialPop = 20; % 5, 10, 20
%L = 50;
thisSet = 1062001;%1066891;%1065301;%1044001;%1064701;%1064501;%1064001;
numItr = 1000%90%90%90;

loadDataOption = 1;
% 1 : load data from full trajectories information
% 2 : load data from sampled trajectories information (*.dat file)

priorConst = 1;%
priorConst2 = 1;%
% pair values [0.1 0.1], [0.1 0.5], [1 1], [1 2]
regStr1 = num2str(priorConst);
regStr2 = num2str(priorConst2);
runLinearInt = 0; % 1: run linear interpolation code, 0: skip it
fileNameContainingDirPath = 'dirNames.txt';
getSysParam;
dT = 10; % 5, 10, 20, 30, 50, 75
ng = 100; % 20, 50, 100, 200
Tstart = 1;
Tused = 100;
%Tend = Tused + Tstart - 1;
Tend = Tused + Tstart;

plotFigs = 0;

timeWholeCodeMPL = zeros(1, numItr);
%%
numParam = Lin*(Lin+1)/2;
allEstEpi = zeros(numItr,numParam);
allEstEpiWithR = zeros(numItr,numParam);
allEstNoMuEpi = zeros(numItr,numParam);
allEstEpiNoE = zeros(numItr,numParam);
allEst = zeros(numItr,Lin);
allNRMSEEpi = zeros(1, numItr);
allNRMSEEpiNoE = zeros(1, numItr);
allNRMSEEpi_Site1n2 = zeros(1, numItr);
allNRMSE = zeros(1, numItr);
allNRMSE_MPL_Site1n2 = zeros(1, numItr);
for thisItr = 1:numItr %2/3, 7, *10*, 190  % 10 is a good example case
    tic

%%  1. load data (in .mat or .dat format)
%--------------------------------------------------------------------------

    thisItr
    mainDir = pwd;    
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

    [dirNameData, dirNameAnalysis] = loadDirNames(fileNameContainingDirPath);
    dirNameRandPermFiles = [dirNameData(1:end-5) chosenSlash 'RandPermFiles' chosenSlash];
    dirNameData = [dirNameData 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
    dirNameAnalysis = [dirNameAnalysis 'Sites' num2str(Lin) chosenSlash selTypeName chosenSlash 'Set' num2str(thisSet) chosenSlash];
    
    
    if(exist(dirNameAnalysis, 'dir') == 0)
        mkdir(dirNameAnalysis)        
    end

%     disp('Loading file...')
%     if(loadDataOption == 1)
%         % load data from full trajectories
%         fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' TgenStr '_' num2str(thisItr) '.mat'];
% 
%         load([dirNameData fileName], 'N', 'D', 'masterStrainList', 'masterTimeFreq', 'perSiteSelction', 'sitesUnderSelection',...
%             'muVal', 'perSiteSelctionEpi_skeleton', 'perSiteSelctionEpi', 'currentNumOfCirStrainsInAllPop', 'freqStatsAllStrains', 'masterSelectionFactor');
% 
% 
%         masterStrainList = masterStrainList(:,1:Lin);
%     elseif(loadDataOption == 2)
%         % load data from sampled trajectories
%         if(Tstart == 1)
%             fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
%         else
%             fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
%         end
%         dataIn = dlmread([dirNameData fileName]);
%     end
disp('Loading file...')
    if(loadDataOption == 1)
        if(recombination == 0)
            % load data from full trajectories
            fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' TgenStr '_' num2str(thisItr) '.mat'];
        elseif(recombination == 1)
            fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_initStr' num2str(numStrainsInInitialPop) '_' TgenStr '_' num2str(thisItr) '.mat'];
        end

        load([dirNameData fileName], 'N', 'D', 'masterStrainList', 'masterTimeFreq', 'perSiteSelction', 'sitesUnderSelection',...
            'muVal', 'perSiteSelctionEpi_skeleton', 'perSiteSelctionEpi', 'currentNumOfCirStrainsInAllPop', 'freqStatsAllStrains', 'masterSelectionFactor');


        masterStrainList = masterStrainList(:,1:Lin);
    elseif(loadDataOption == 2)
        if(recombination == 0)
            % load data from sampled trajectories
            if(Tstart == 1)
                fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
            else
                fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
            end
        elseif(recombination == 1)
            % load data from sampled trajectories
            if(Tstart == 1)
                fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
            else
                fileName = ['WFsimEpi_Nmu' NmuInStr '_Nr' NrInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_Dele' num2str(DeleAll) '_selVal' num2str(selVal) '_delSelVal' num2str(delSelVal) '_' TgenStr '_' num2str(thisItr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
            end
        end
        dataIn = dlmread([dirNameData fileName]);
    end
    
    warning off
%% 2. calculates one and two point frequencies
%--------------------------------------------------------------------------
    if(loadDataOption == 1)
        % Reconstruct MSA from raw data files, perform finite sampling by 
        % selecting ng individuals from the MSA and find sampled 1 and 2 point
        % frequencies
        %numGen = (Tend - Tstart + 1)/dT; % number of generation in the whole msa
        numGen = (Tend - Tstart)/dT + 1; % number of generation in the whole msa

        % later can make a switch to chose less than N samples too to simulate what
        % happens in FLU samples 
        numSamplesPerGen = N; 
        numSamplesPerGenSelected = ng;

        samplingTimePoints = Tstart:dT:Tend;
        numSamplingPoints = length(samplingTimePoints);
        q = -1*ones(numSamplingPoints, Lin);
        q11 = -1*ones(Lin, Lin, numSamplingPoints);
                
        numUnique2PtPairs = (Lin*(Lin-1)/2);
        qTwoPt_jk_Identity = zeros(numUnique2PtPairs, 2);
        countTemp1 = 1;
        for l2 = 1:Lin-1
            for l3 = l2+1:Lin
                qTwoPt_jk_Identity(countTemp1,:) = [l2 l3];
                countTemp1 = countTemp1 + 1;
            end
        end
        
        sumAijWithLink = zeros(Lin,Lin);
        sumAijWithLink2 = zeros(Lin,Lin);
        sumEWithLink = zeros(Lin, numUnique2PtPairs);
        sumFWithLink = zeros(numUnique2PtPairs, numUnique2PtPairs);
        qExt = zeros(numSamplingPoints, Lin + numUnique2PtPairs);
%         CmatAll = zeros(5, 5, numSamplingPoints);
        HmatAll = zeros(numParam, numParam, numSamplingPoints);
        
        saved4pt = -1*ones(1, numSamplingPoints);
        %randSelctIndAll = zeros(numSamplingPoints, numSamplesPerGenSelected);
        % this section calculates q_t for each site, for the synthetic protein, 
        % 2 allele per site case this is just the sum of occurances of ones
        fprintf('Calculate 1 and 2 point probabilities...')
        for t = 1:numSamplingPoints
            
            thisSamplingTimePoint = samplingTimePoints(t);
            strainsThisTimePoint = masterTimeFreq{thisSamplingTimePoint}(:,1);
            freqStrainThisTimePoint = masterTimeFreq{thisSamplingTimePoint}(:,2);
            thisMSATemp = -1*ones(numSamplesPerGen, Lin);
            q11Temp = -1*ones(Lin, Lin);
            temp200 = zeros(Lin, Lin);
            temp210 = zeros(Lin, Lin);
            temp300 = zeros(Lin, numUnique2PtPairs);
            temp400 = zeros(numUnique2PtPairs, numUnique2PtPairs);
            count1 = 0;
            for k = 1:length(strainsThisTimePoint)
                count1 = count1 + 1;
                thisMSATemp(count1:count1+freqStrainThisTimePoint(k)-1,:) = repmat(masterStrainList(strainsThisTimePoint(k),:), freqStrainThisTimePoint(k), 1);
                count1 = count1 + freqStrainThisTimePoint(k) - 1;
            end
            fileNameRandPerm = ['RandPerm_N' num2str(N) '_t' num2str(t)];
            load([dirNameRandPermFiles fileNameRandPerm],'randPermN')
            temp4051 = randPermN;
            randSelctInd = temp4051(1:numSamplesPerGenSelected);
            %randSelctIndAll(t, :) = randSelctInd;

            thisMSA = thisMSATemp(randSelctInd,:);
            q(t,:) = sum(thisMSA)./numSamplesPerGenSelected; % normalize frequencies;
            thisMSALogical = logical(thisMSA);
            tempDiagEntries = q(t,:).*(1 - q(t,:));
            allSiteInd = 1:Lin;
            qExtThis = zeros(1,Lin + numUnique2PtPairs);
            qExtThis(1:Lin) = q(t,:);
            count = Lin;
            for l = 1:Lin
                % multiply the lth column of MSA with all MSA to find q11 for (l,:)
                % pairs of cov matrix
                tempq11 = (repmat(thisMSALogical(:,l),1, Lin).* thisMSALogical);
                
                % sum for the lth row of the covariance matrix
                %q11(l,:,t) = sum(tempq11)./numSamplesPerGenSelected;
                
                q11Temp(l,:) = sum(tempq11)./numSamplesPerGenSelected;
                q11(l,:,t) = q11Temp(l,:);
                thisEntryLength = Lin - l;
                qExtThis(count+1:count + thisEntryLength) = q11Temp(l,l+1:end);
                count = count + thisEntryLength;
                
                temp100 = q(t,:).*q(t,l);
                temp200(l,:) = q11Temp(l, :) - temp100;
                                
                % calculate 3-points and calculate q111-q1*q11                
                q111Temp_partial = zeros(1, numUnique2PtPairs);
                q11_jkTemp_partial = zeros(1, numUnique2PtPairs);
                for ll = 1:numUnique2PtPairs
                    jthSite = qTwoPt_jk_Identity(ll,1);
                    kthSite = qTwoPt_jk_Identity(ll,2);
                    if(l <= Lin)
                        q111Temp_partial(ll) = sum(thisMSALogical(:,l).*thisMSALogical(:,jthSite).*thisMSALogical(:,kthSite))/numSamplesPerGenSelected;
                        q11_jkTemp_partial(ll) = sum(thisMSALogical(:,jthSite).*thisMSALogical(:,kthSite))/numSamplesPerGenSelected;
                    end
                end

                temp300(l,:) = q111Temp_partial - q11_jkTemp_partial.*q(t,l);

            end
            qExtThis(1:Lin) = q(t,:)';
            qExtThis(Lin+1:end) = q11_jkTemp_partial';
            qExt(t,:) = qExtThis;
            % calc 4 points
            for l2 = 1:numUnique2PtPairs
                hthSite = qTwoPt_jk_Identity(l2,1);
                ithSite = qTwoPt_jk_Identity(l2,2);
                q11_hiTemp = sum(thisMSALogical(:,hthSite).*thisMSALogical(:,ithSite))/numSamplesPerGenSelected;
                q11_11Temp_partial = zeros(1, numUnique2PtPairs);
                q11_jkTemp_partial = zeros(1, numUnique2PtPairs);
                for ll = 1:numUnique2PtPairs
                    jthSite = qTwoPt_jk_Identity(ll,1);
                    kthSite = qTwoPt_jk_Identity(ll,2);

                    q11_11Temp_partial(ll) = sum(thisMSALogical(:,hthSite).*thisMSALogical(:,ithSite).*thisMSALogical(:,jthSite).*thisMSALogical(:,kthSite))/numSamplesPerGenSelected;
                    q11_jkTemp_partial(ll) = sum(thisMSALogical(:,jthSite).*thisMSALogical(:,kthSite))/numSamplesPerGenSelected;
                    
                    if(hthSite == 1 && ithSite == 2 && jthSite == 3 && kthSite == 4)
                        %[t l ll]
                        saved4pt(t) = q11_11Temp_partial(ll);
                        %pause
                    end
                end
                temp400(l2,:) = q11_11Temp_partial - q11_hiTemp*q11_jkTemp_partial;
                
            end
            
            if(t == 1)
                qExt_tStart = qExtThis;
            elseif(t == numSamplingPoints)
                qExt_tLast = qExtThis;
            end
            
            covMtxThisTime = ~(eye(Lin)).*temp200 + diag(tempDiagEntries);
            
            if(t < numSamplingPoints)
                sumAijWithLink = sumAijWithLink + covMtxThisTime;

                sumAijWithLink2 = sumAijWithLink2 + temp200;

                sumEWithLink = sumEWithLink + temp300;
                sumFWithLink = sumFWithLink + temp400;
                
                HmatAll(:,:,t) = [covMtxThisTime temp300;
                                  temp300' temp400];
            end
        end
        
        perSiteSelctionAllEpiTerms = [];
        for l = 1:Lin
            perSiteSelctionAllEpiTerms = [perSiteSelctionAllEpiTerms perSiteSelctionEpi(l,l+1:end)];
        end

        disp('done')

        clear qijAtTimeTk;

        clear masterStrainFreq;
    elseif(loadDataOption == 2)
        AllMSATimeVec = dataIn(:,1) + 1;
        AllMSAFreqIn = dataIn(:,2);
        AllMSAIn = dataIn(:,3:end);
        samplingTimePoints = unique(dataIn(:,1) + 1)';
        numSamplingPoints = length(samplingTimePoints);
        numSamplesPerGen = Nin; 
        numSamplesPerGenSelected = ng;
        numGen = length(samplingTimePoints);
        q = -1*ones(numSamplingPoints, Lin);
        %q11 = -1*ones(Lin, Lin, numSamplingPoints);
        sumAijWithLink = zeros(Lin,Lin);
        % this section calculates q_t for each site, for the synthetic protein, 
        % 2 allele per site case this is just the sum of occurances of ones
        fprintf('Calculate 1 and 2 point probabilities...')
        for t = 1:numSamplingPoints

            thisSamplingTimePoint = samplingTimePoints(t);
            thisTimePointSelcRows = AllMSATimeVec == thisSamplingTimePoint;
            strainsThisTimePoint = AllMSAIn(thisTimePointSelcRows,:);
            freqStrainThisTimePoint = AllMSAFreqIn(thisTimePointSelcRows);
            thisMSA = -1*ones(numSamplesPerGenSelected, Lin);
            q11Temp = -1*ones(Lin, Lin);
            temp200 = zeros(Lin,Lin);
            
            count1 = 0;
            for k = 1:size(strainsThisTimePoint,1)
                count1 = count1 + 1;
                thisMSA(count1:count1+freqStrainThisTimePoint(k)-1,:) = repmat(strainsThisTimePoint(k,:), freqStrainThisTimePoint(k), 1);
                count1 = count1 + freqStrainThisTimePoint(k) - 1;
            end

            q(t,:) = sum(thisMSA)./numSamplesPerGenSelected; % normalize frequencies;
            thisMSALogical = logical(thisMSA);
            tempDiagEntries = q(t,:).*(1 - q(t,:));
            
            for l = 1:Lin
                % multiply the lth column of MSA with all MSA to find q11 for (l,:)
                % pairs of cov matrix
                tempq11 = (repmat(thisMSALogical(:,l),1, Lin).* thisMSALogical);
                % sum for the lth row of the covariance matrix
                %q11(l,:,t) = sum(tempq11)./numSamplesPerGenSelected;
                q11Temp(l,:) = sum(tempq11)./numSamplesPerGenSelected;
            
                temp100 = q(t,:).*q(t,l);
                temp200(l,:) = q11Temp(l, :) - temp100;
            end            
            covMtxThisTime = ~(eye(Lin)).*temp200 + diag(tempDiagEntries);
            sumAijWithLink = sumAijWithLink + covMtxThisTime;
        end
        disp('done')
    end
    
    % plot tranjectory at each site
    if(plotFigs)
        for l = 1:Lin
           figure(l)
           subplot(2,1,1)
           plot(q(:,l), 'b.-')
           xlabel('time points')
           ylabel('Frequency Count')
           title(['Site number: ' num2str(l)])
        end
    end
    
%% 3. Find estimates using MPL
%     vectors with mu
    
    fprintf('Calculating selection coefficients estimates...')
    model = 'linkDiff';

    startTime = 1;
    stopTime = numGen;
    
    tempTimeSumCovProcessig = toc;

    % 6.1.2.1 calc estimates of selection and errors
    sumAijUnLink = diag(diag(sumAijWithLink));

    v_EstOutLink = dT*muVal*(sum((ones(Lin,stopTime-startTime+1-1) - 2*q(startTime:stopTime-1,:)'),2));
    
    temp40 = zeros(numUnique2PtPairs, 1);
    for t = startTime:stopTime - 1
        count3 = 0;
        for l2 = 1:Lin-1
            temp40(count3 + 1:count3 + Lin - l2) = temp40(count3 + 1:count3 + Lin - l2) + repmat(q(t,l2), Lin - l2, 1) + q(t, l2+1:end)';
            count3 = count3 + Lin - l2;
        end
    end
    
    temp60 = dT*muVal*(temp40-4*sum(qExt(startTime:stopTime-1,Lin+1:end))');

    %temp60 = dT*muVal*sum(temp40-4*qExt(startTime:stopTime-1,Lin+1:end)',2);
    vExt_EstOutLink = [v_EstOutLink; temp60];
    
    Reg1Mat = priorConst*eye(Lin);
    Reg2Mat = diag([priorConst*ones(1,Lin) priorConst2*ones(1, numUnique2PtPairs)]);
    
    sigmaEstOutLink = (sumAijWithLink*dT + Reg1Mat)\(q(stopTime,:) - q(startTime,:) - v_EstOutLink')';
    sigmaEstOutLinkNoMu = (sumAijWithLink*dT + Reg1Mat)\(q(stopTime,:) - q(startTime,:))';
    sigmaEstOutUnLink = (sumAijUnLink*dT + Reg1Mat)\(q(stopTime,:) - q(startTime,:) - v_EstOutLink')';
    sigmaEstOutUnLinkNoMu = (sumAijUnLink*dT + Reg1Mat)\(q(stopTime,:) - q(startTime,:))';

    Hmat = [sumAijWithLink sumEWithLink;
            sumEWithLink'  sumFWithLink];
        
    if(useStrongReg == 1)
        [selcSitesEpiAcc, ~,~,~] = findIndColsOfHmat(Hmat*dT);
        selcSitesEpiInAcc = ~selcSitesEpiAcc;
        selcSitesEpiInAccEpiTerms = selcSitesEpiInAcc(Lin+1:end);
        priorConst2Vec = priorConst2*ones(1, numUnique2PtPairs);
        priorConst2Vec(selcSitesEpiInAccEpiTerms) = priorConst2*10;
    end
    
    
    
    sigmaEstOutLinkEpi = (Hmat*dT + Reg2Mat)\(qExt_tLast - qExt_tStart - vExt_EstOutLink')';
    sigmaEstOutLinkNoMuEpi = inv(Hmat*dT + Reg2Mat)*(qExt_tLast - qExt_tStart)';

    HmatNoE = [sumAijWithLink sumEWithLink*0;
            0*sumEWithLink'  sumFWithLink];
    sigmaEstOutLinkEpiNoE = (HmatNoE*dT + Reg2Mat)\(qExt_tLast - qExt_tStart - vExt_EstOutLink')';
    sigmaEstOutLinkNoMuEpiNoE = inv(HmatNoE*dT + Reg2Mat)*(qExt_tLast - qExt_tStart)';

    
    if(recombination == 1)
        recTwoPointsTemp = zeros(Lin*(Lin - 1)/2,1); % only epi trms will have this
        count1 = 1;
        for rr = 1:Lin-1
            mulVec = 1:Lin-rr; % mulVec is magnitude of (i-j)
            numEntriesThis = length(mulVec);
            temprr2 = mulVec.*sumAijWithLink(rr, rr+1:Lin);
            recTwoPointsTemp(count1:count1+numEntriesThis-1) = temprr2';
            count1 = count1 + numEntriesThis;
        end
        temp33 = recTwoPointsTemp;
        uExt_EstOutLinkWithR = -dT*recVal*[zeros(Lin, 1); temp33]; % (i-j) is mulVec
        numLinkEpiWithR = (qExt_tLast - qExt_tStart - vExt_EstOutLink' - uExt_EstOutLinkWithR')';

        sigmaEstOutLinkEpiWithR = (Hmat*dT + Reg2Mat)\numLinkEpiWithR;
    else
        sigmaEstOutLinkEpiWithR = zeros(Lin*(Lin+1)/2, 1);
    end
% figure
% plot(1:15, diag(inv(Hmat*dT + Reg2Mat))', 'bo-')
% hold on
% %plot(1:15, diag(Hmat*dT)', 'r.-')
% plot(1:15, diag(Hmat*dT + Reg2Mat)'/10, 'r.-')
% pause
    % likelihood stuff needs to be checked to make it consistent with
    % different reg for epi terms
%------------- likelihood, vec EPi length 15, MPl, SL length 5-----------
    regConst2 = 0;
    logLikelihoodEpi = 0;
    logLikelihoodEpi0 = 0;
    logLikelihoodMPL = 0;
    logLikelihoodMPL0 = 0;
    logLikelihoodSL = 0;
    logLikelihood0 = 0;
    timePointsUsedNoEpi = 0;
    timePointsUsedEpi = 0;
    numParamEstNoEpi = 0;
    numParamEstEpi = 0;
    
    xCount = 0;
    yCount = 0;
    logLikeliTempSL_term1 = zeros(1,numSamplingPoints);
    logLikeliTempSL_term2 = zeros(1,numSamplingPoints);
    logLikeliTempSL_term3 = zeros(1,numSamplingPoints);
    logLikeliTempSL_term4 = zeros(1,numSamplingPoints);
    logLikeliTempMPL_term1 = zeros(1,numSamplingPoints);
    logLikeliTempMPL_term2 = zeros(1,numSamplingPoints);
    logLikeliTempMPL_term3 = zeros(1,numSamplingPoints);
    logLikeliTempMPL_term4 = zeros(1,numSamplingPoints);
    logLikeliTempMPL0_term1 = zeros(1,numSamplingPoints);
    logLikeliTempMPL0_term2 = zeros(1,numSamplingPoints);
    logLikeliTempMPL0_term3 = zeros(1,numSamplingPoints);
    logLikeliTempMPL0_term4 = zeros(1,numSamplingPoints);
    logLikeliEpiTemp_term1 = zeros(1,numSamplingPoints);
    logLikeliEpiTemp_term2 = zeros(1,numSamplingPoints);
    logLikeliEpiTemp_term3 = zeros(1,numSamplingPoints);
    logLikeliEpiTemp_term4 = zeros(1,numSamplingPoints);
    logLikeliEpiTemp0_term1 = zeros(1,numSamplingPoints);
    logLikeliEpiTemp0_term2 = zeros(1,numSamplingPoints);
    logLikeliEpiTemp0_term3 = zeros(1,numSamplingPoints);
    logLikeliEpiTemp0_term4 = zeros(1,numSamplingPoints);
    
    for t = 1:numSamplingPoints-1
        % making the V vetor for the Epi case
        temp40 = zeros(numUnique2PtPairs, 1);

        count3 = 0;
        for l2 = 1:Lin-1
            temp40(count3 + 1:count3 + Lin - l2) = temp40(count3 + 1:count3 + Lin - l2) + repmat(q(t,l2), Lin - l2, 1) + q(t, l2+1:end)';
            count3 = count3 + Lin - l2;
        end

        vEpi = [ (1 - 2*q(t,:)');  (temp40-4*qExt(t,Lin+1:end)')];
        thisHmat = HmatAll(:,:,t);
        [RTemp, jbThisHmat] = rref(thisHmat);
        thisCmat = thisHmat(1:Lin,1:Lin);
        [RTemp, jbThisCmat] = rref(thisCmat);
        if(~isempty(jbThisCmat))
            xCount = xCount + 1;
            timePointsUsedNoEpi = timePointsUsedNoEpi + 1;
            numParamEstNoEpi = numParamEstNoEpi + length(jbThisCmat);
            thisCmatUsable = thisCmat(jbThisCmat,jbThisCmat);

            % SL SC
            thisCmatSL = diag(diag(thisCmatUsable));
            denSL = det(thisCmatSL + regConst2*eye(length(jbThisCmat)));
            vecSL = q(t+1,jbThisCmat)' - q(t,jbThisCmat)' - dT*(thisCmatSL + regConst2*eye(length(jbThisCmat)))*sigmaEstOutUnLink(jbThisCmat) - dT*muVal*(1 - 2*q(t,jbThisCmat)');
            %logLikeliTempSL = -Nin/(4*dT)*vecSL'*(thisCmatSL + priorConst*eye(length(jbThisCmat)))^-1*vecSL - 1/2*log(denSL) + length(jbThisCmat)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
            logLikeliTempSL_term1(xCount) = -Nin/(4*dT)*vecSL'*(thisCmatSL + 0*priorConst*eye(length(jbThisCmat)))^-1*vecSL;
            logLikeliTempSL_term2(xCount) = - 1/2*log(denSL);
            logLikeliTempSL_term3(xCount) = length(jbThisCmat)/2*log(Nin/(4*pi*dT));
            logLikeliTempSL_term4(xCount) = (log(sum(1/Nin*length(jbThisCmat))));%(log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
            %logLikeliTempSL = logLikeliTempSL_term1(xCount) + logLikeliTempSL_term2(xCount) + logLikeliTempSL_term4(xCount);
            logLikeliTempSL = logLikeliTempSL_term1(xCount) + logLikeliTempSL_term2(xCount) + logLikeliTempSL_term3(xCount) + logLikeliTempSL_term4(xCount);
            logLikelihoodSL = logLikelihoodSL + logLikeliTempSL;

            % MPL SC
            den = det(thisCmatUsable + regConst2*eye(length(jbThisCmat)));
            vec = q(t+1,jbThisCmat)' - q(t,jbThisCmat)' - dT*(thisCmatUsable + regConst2*eye(length(jbThisCmat)))*sigmaEstOutLink(jbThisCmat) - dT*muVal*(1 - 2*q(t,jbThisCmat)');
            logLikeliTempMPL_term1(xCount) = -Nin/(4*dT)*vec'*(thisCmatUsable + 0*priorConst*eye(length(jbThisCmat)))^-1*vec;
            logLikeliTempMPL_term2(xCount) = - 1/2*log(den);
            logLikeliTempMPL_term3(xCount) = length(jbThisCmat)/2*log(Nin/(4*pi*dT));
            logLikeliTempMPL_term4(xCount) = (log(sum(1/Nin*length(jbThisCmat))));%(log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
            %logLikeliTempMPL = logLikeliTempMPL_term1(xCount) + logLikeliTempMPL_term2(xCount) + logLikeliTempMPL_term4(xCount);
            logLikeliTempMPL = logLikeliTempMPL_term1(xCount) + logLikeliTempMPL_term2(xCount) + logLikeliTempMPL_term3(xCount) + logLikeliTempMPL_term4(xCount);
            %logLikeliTempMPL = -Nin/(4*dT)*vec'*(thisCmatUsable + priorConst*eye(length(jbThisCmat)))^-1*vec - 1/2*log(den) + length(jbThisCmat)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
            logLikelihoodMPL = logLikelihoodMPL + logLikeliTempMPL;
            
            
            % MPL with 0 SC
            den0 = det(thisCmatUsable + regConst2*eye(length(jbThisCmat)));
            vec0 = q(t+1,jbThisCmat)' - q(t,jbThisCmat)' - dT*(thisCmatUsable + regConst2*eye(length(jbThisCmat)))*sigmaEstOutLink(jbThisCmat)*0 - dT*muVal*(1 - 2*q(t,jbThisCmat)');
            logLikeliTempMPL0_term1(xCount) = -Nin/(4*dT)*vec0'*(thisCmatUsable + 0*priorConst*eye(length(jbThisCmat)))^-1*vec0;
            logLikeliTempMPL0_term2(xCount) = - 1/2*log(den0);
            logLikeliTempMPL0_term3(xCount) = length(jbThisCmat)/2*log(Nin/(4*pi*dT));
            logLikeliTempMPL0_term4(xCount) = (log(sum(1/Nin*length(jbThisCmat))));%(log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
            %logLikeliTempMPL = logLikeliTempMPL_term1(xCount) + logLikeliTempMPL_term2(xCount) + logLikeliTempMPL_term4(xCount);
            logLikeliTempMPL0 = logLikeliTempMPL0_term1(xCount) + logLikeliTempMPL0_term2(xCount) + logLikeliTempMPL0_term3(xCount) + logLikeliTempMPL0_term4(xCount);
            %logLikeliTempMPL = -Nin/(4*dT)*vec'*(thisCmatUsable + priorConst*eye(length(jbThisCmat)))^-1*vec - 1/2*log(den) + length(jbThisCmat)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
            logLikelihoodMPL0 = logLikelihoodMPL0 + logLikeliTempMPL0;
        end

        if(~isempty(jbThisHmat))
            yCount = yCount + 1;
            timePointsUsedEpi = timePointsUsedEpi + 1;
            numParamEstEpi = numParamEstEpi + length(jbThisHmat);
            thisHmatUsable = thisHmat(jbThisHmat,jbThisHmat);
            %MPLE SC
            denEpi = det(thisHmatUsable + regConst2*eye(length(jbThisHmat)));
            vecEpi = qExt(t+1,jbThisHmat)' - qExt(t,jbThisHmat)' - dT*thisHmatUsable*sigmaEstOutLinkEpi(jbThisHmat) - dT*muVal*vEpi(jbThisHmat);
            logLikeliEpiTemp_term1(yCount) = -Nin/(4*dT)*vecEpi'*(thisHmatUsable + 0*priorConst*eye(length(jbThisHmat)))^-1*vecEpi;
            logLikeliEpiTemp_term2(yCount) = - 1/2*log(denEpi);
            logLikeliEpiTemp_term3(yCount) = length(jbThisHmat)/2*log(Nin/(4*pi*dT));
            logLikeliEpiTemp_term4(yCount) = (log(sum(1/Nin*length(jbThisHmat))));%(log(sum(abs(qExt(t+1,jbThisHmat)' - qExt(t,jbThisHmat)'))));
            %logLikeliEpiTemp = logLikeliEpiTemp_term1(yCount) + logLikeliEpiTemp_term2(yCount) + logLikeliEpiTemp_term4(yCount);
            logLikeliEpiTemp = logLikeliEpiTemp_term1(yCount) + logLikeliEpiTemp_term2(yCount) + logLikeliEpiTemp_term3(yCount) + logLikeliEpiTemp_term4(yCount);
            %logLikeliEpiTemp = -Nin/(4*dT)*vecEpi'*(thisHmatUsable + priorConst*eye(length(jbThisHmat)))^-1*vecEpi - 1/2*log(denEpi) + (length(jbThisHmat)*(length(jbThisHmat) + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(qExt(t+1,jbThisHmat)' - qExt(t,jbThisHmat)'))));
            logLikelihoodEpi = logLikelihoodEpi + logLikeliEpiTemp;
            
            % MPLE 0 SC
            denEpi0 = det(thisHmatUsable + regConst2*eye(length(jbThisHmat)));
            vecEpi0 = qExt(t+1,jbThisHmat)' - qExt(t,jbThisHmat)' - dT*thisHmatUsable*sigmaEstOutLinkEpi(jbThisHmat)*0 - dT*muVal*vEpi(jbThisHmat);
            logLikeliEpiTemp0_term1(yCount) = -Nin/(4*dT)*vecEpi0'*(thisHmatUsable + 0*priorConst*eye(length(jbThisHmat)))^-1*vecEpi0;
            logLikeliEpiTemp0_term2(yCount) = - 1/2*log(denEpi0);
            logLikeliEpiTemp0_term3(yCount) = length(jbThisHmat)/2*log(Nin/(4*pi*dT));
            logLikeliEpiTemp0_term4(yCount) = (log(sum(1/Nin*length(jbThisHmat))));%(log(sum(abs(qExt(t+1,jbThisHmat)' - qExt(t,jbThisHmat)'))));
            %logLikeliEpiTemp = logLikeliEpiTemp_term1(yCount) + logLikeliEpiTemp_term2(yCount) + logLikeliEpiTemp_term4(yCount);
            logLikeliEpiTemp0 = logLikeliEpiTemp0_term1(yCount) + logLikeliEpiTemp0_term2(yCount) + logLikeliEpiTemp0_term3(yCount) + logLikeliEpiTemp0_term4(yCount);
            %logLikeliEpiTemp = -Nin/(4*dT)*vecEpi'*(thisHmatUsable + priorConst*eye(length(jbThisHmat)))^-1*vecEpi - 1/2*log(denEpi) + (length(jbThisHmat)*(length(jbThisHmat) + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(qExt(t+1,jbThisHmat)' - qExt(t,jbThisHmat)'))));
            logLikelihoodEpi0 = logLikelihoodEpi0 + logLikeliEpiTemp0;
        end
        
        [logLikelihoodSL logLikelihoodMPL logLikelihoodMPL0 logLikelihoodEpi logLikelihoodEpi0]
        
        temp4040 = -2*[logLikelihoodSL logLikelihoodMPL logLikelihoodMPL0 logLikelihoodEpi logLikelihoodEpi0] + [5 5 5 15 15].*log([ng*timePointsUsedEpi]);
        %temp4040 = -2*[logLikelihoodSL logLikelihoodMPL logLikelihoodEpi] + timePointsUsedNoEpi*[numParamEstNoEpi numParamEstNoEpi numParamEstEpi].*log([numParamEstNoEpi*timePointsUsedNoEpi numParamEstNoEpi*timePointsUsedNoEpi numParamEstEpi*timePointsUsedEpi]);
        


% %------------- likelihood, vec EPi length 15, MPl, SL length 5-----------
%     regConst2 = 0;
%     logLikelihoodEpi = 0;
%     logLikelihoodMPL = 0;
%     logLikelihoodSL = 0;
%     logLikelihood0 = 0;
%     timePointsUsedNoEpi = 0;
%     timePointsUsedEpi = 0;
%     numParamEstNoEpi = 0;
%     numParamEstEpi = 0;
%     xCount = 0;
%     yCount = 0;
%     for t = 1:numSamplingPoints-1
%         % making the V vetor for the Epi case
%         temp40 = zeros(numUnique2PtPairs, 1);
% 
%         count3 = 0;
%         for l2 = 1:Lin-1
%             temp40(count3 + 1:count3 + Lin - l2) = temp40(count3 + 1:count3 + Lin - l2) + repmat(q(t,l2), Lin - l2, 1) + q(t, l2+1:end)';
%             count3 = count3 + Lin - l2;
%         end
% 
%         vEpi = [ (1 - 2*q(t,:)');  (temp40-4*qExt(t,Lin+1:end)')];
%         thisHmat = HmatAll(:,:,t);
%         [RTemp, jbThisHmat] = rref(thisHmat);
%         thisCmat = thisHmat(1:Lin,1:Lin);
%         [RTemp, jbThisCmat] = rref(thisCmat);
%         if(~isempty(jbThisCmat))
%             xCount = xCount + 1;
%             timePointsUsedNoEpi = timePointsUsedNoEpi + 1;
%             numParamEstNoEpi = numParamEstNoEpi + length(jbThisCmat);
%             thisCmatUsable = thisCmat(jbThisCmat,jbThisCmat);
% 
%             % SL SC
%             thisCmatSL = diag(diag(thisCmatUsable));
%             denSL = det(thisCmatSL + regConst2*eye(length(jbThisCmat)));
%             vecSL = q(t+1,jbThisCmat)' - q(t,jbThisCmat)' - dT*(thisCmatSL + regConst2*eye(length(jbThisCmat)))*sigmaEstOutUnLink(jbThisCmat) - dT*muVal*(1 - 2*q(t,jbThisCmat)');
%             %logLikeliTempSL = -Nin/(4*dT)*vecSL'*(thisCmatSL + priorConst*eye(length(jbThisCmat)))^-1*vecSL - 1/2*log(denSL) + length(jbThisCmat)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
%             logLikeliTempSL_term1(xCount) = -Nin/(4*dT)*vecSL'*(thisCmatSL + 0*priorConst*eye(length(jbThisCmat)))^-1*vecSL;
%             logLikeliTempSL_term2(xCount) = - 1/2*log(denSL);
%             logLikeliTempSL_term3(xCount) = length(jbThisCmat)/2*log(Nin/(4*pi*dT));
%             logLikeliTempSL_term4(xCount) = (log(sum(1/Nin*length(jbThisCmat))));%(log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
%             %logLikeliTempSL = logLikeliTempSL_term1(xCount) + logLikeliTempSL_term2(xCount) + logLikeliTempSL_term4(xCount);
%             logLikeliTempSL = logLikeliTempSL_term1(xCount) + logLikeliTempSL_term2(xCount) + logLikeliTempSL_term3(xCount) + logLikeliTempSL_term4(xCount);
%             logLikelihoodSL = logLikelihoodSL + logLikeliTempSL;
% 
%             % MPL SC
%             den = det(thisCmatUsable + regConst2*eye(length(jbThisCmat)));
%             vec = q(t+1,jbThisCmat)' - q(t,jbThisCmat)' - dT*(thisCmatUsable + regConst2*eye(length(jbThisCmat)))*sigmaEstOutLink(jbThisCmat) - dT*muVal*(1 - 2*q(t,jbThisCmat)');
%             logLikeliTempMPL_term1(xCount) = -Nin/(4*dT)*vec'*(thisCmatUsable + 0*priorConst*eye(length(jbThisCmat)))^-1*vec;
%             logLikeliTempMPL_term2(xCount) = - 1/2*log(den);
%             logLikeliTempMPL_term3(xCount) = length(jbThisCmat)/2*log(Nin/(4*pi*dT));
%             logLikeliTempMPL_term4(xCount) = (log(sum(1/Nin*length(jbThisCmat))));%(log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
%             %logLikeliTempMPL = logLikeliTempMPL_term1(xCount) + logLikeliTempMPL_term2(xCount) + logLikeliTempMPL_term4(xCount);
%             logLikeliTempMPL = logLikeliTempMPL_term1(xCount) + logLikeliTempMPL_term2(xCount) + logLikeliTempMPL_term3(xCount) + logLikeliTempMPL_term4(xCount);
%             %logLikeliTempMPL = -Nin/(4*dT)*vec'*(thisCmatUsable + priorConst*eye(length(jbThisCmat)))^-1*vec - 1/2*log(den) + length(jbThisCmat)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(q(t+1,jbThisCmat)' - q(t,jbThisCmat)'))));
%             logLikelihoodMPL = logLikelihoodMPL + logLikeliTempMPL;
%         end
% 
%         if(~isempty(jbThisHmat))
%             yCount = yCount + 1;
%             timePointsUsedEpi = timePointsUsedEpi + 1;
%             numParamEstEpi = numParamEstEpi + length(jbThisHmat);
%             thisHmatUsable = thisHmat(jbThisHmat,jbThisHmat);
%             %MPLE SC
%             denEpi = det(thisHmatUsable + regConst2*eye(length(jbThisHmat)));
%             vecEpi = qExt(t+1,jbThisHmat)' - qExt(t,jbThisHmat)' - dT*thisHmatUsable*sigmaEstOutLinkEpi(jbThisHmat)*3 - dT*muVal*vEpi(jbThisHmat);
%             logLikeliEpiTemp_term1(yCount) = -Nin/(4*dT)*vecEpi'*(thisHmatUsable + 0*priorConst*eye(length(jbThisHmat)))^-1*vecEpi;
%             logLikeliEpiTemp_term2(yCount) = - 1/2*log(denEpi);
%             logLikeliEpiTemp_term3(yCount) = length(jbThisHmat)/2*log(Nin/(4*pi*dT));
%             logLikeliEpiTemp_term4(yCount) = (log(sum(1/Nin*length(jbThisCmat))));%(log(sum(abs(qExt(t+1,jbThisHmat)' - qExt(t,jbThisHmat)'))));
%             %logLikeliEpiTemp = logLikeliEpiTemp_term1(yCount) + logLikeliEpiTemp_term2(yCount) + logLikeliEpiTemp_term4(yCount);
%             logLikeliEpiTemp = logLikeliEpiTemp_term1(yCount) + logLikeliEpiTemp_term2(yCount) + logLikeliEpiTemp_term3(yCount) + logLikeliEpiTemp_term4(yCount);
%             %logLikeliEpiTemp = -Nin/(4*dT)*vecEpi'*(thisHmatUsable + priorConst*eye(length(jbThisHmat)))^-1*vecEpi - 1/2*log(denEpi) + (length(jbThisHmat)*(length(jbThisHmat) + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(qExt(t+1,jbThisHmat)' - qExt(t,jbThisHmat)'))));
%             logLikelihoodEpi = logLikelihoodEpi + logLikeliEpiTemp;
%         end
%         
%         [logLikelihoodSL logLikelihoodMPL logLikelihoodEpi]
%         
%         temp4040 = -2*[logLikelihoodSL logLikelihoodMPL logLikelihoodEpi] + [5 5 15].*log([numParamEstNoEpi*timePointsUsedNoEpi numParamEstNoEpi*timePointsUsedNoEpi numParamEstEpi*timePointsUsedEpi]);
%         %temp4040 = -2*[logLikelihoodSL logLikelihoodMPL logLikelihoodEpi] + timePointsUsedNoEpi*[numParamEstNoEpi numParamEstNoEpi numParamEstEpi].*log([numParamEstNoEpi*timePointsUsedNoEpi numParamEstNoEpi*timePointsUsedNoEpi numParamEstEpi*timePointsUsedEpi]);
%         
% %--------------------------------------------------------------------------
        
        
        
        %%%% very old liklehood *****wrong*****
%             
%             % MPL SC
%             den = det(thisCmat + regConst2*eye(Lin));
%             vec = q(t+1,:)' - q(t,:)' - dT*(thisCmat + regConst2*eye(Lin))*sigmaEstOutLink - dT*muVal*(1 - 2*q(t,:)');
%             logLikeliTemp = -Nin/(4*dT)*vec'*(thisCmat + regConst2*eye(Lin))^-1*vec - 1/2*log(den) + Lin/2*log(Nin/(4*pi*dT)) + (log(sum(abs(q(t+1,:)' - q(t,:)'))))
% %             den = det(thisHmat + regConst2*eye(15));
% %             vec = qExt(t+1,:)' - qExt(t,:)' - dT*thisHmat*[sigmaEstOutLink; zeros(10,1)] - dT*muVal*vEpi;
% %             logLikeliTemp = -Nin/(4*dT)*vec'*(thisHmat + regConst2*eye(15))^-1*vec - 1/2*log(den) + (Lin*(Lin + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(qExt(t+1,:)' - qExt(t,:)'))));
% %             thisHmatMPL = [thisCmat zeros(5,10); zeros(10,15)];
% %             den = det(thisHmatMPL + regConst2*eye(15));
% %             vec = [qExt(t+1,1:Lin) zeros(1,10)]' - [qExt(t,1:Lin) zeros(1,10)]' - dT*thisHmatMPL*[sigmaEstOutUnLink; zeros(10,1)] - dT*muVal*[vEpi(1:Lin); zeros(10,1)];
% %             logLikeliTemp = -Nin/(4*dT)*vec'*(thisHmatMPL + regConst2*eye(15))^-1*vec - 1/2*log(den) + (Lin*(Lin + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs([qExt(t+1,1:Lin) zeros(1,10)]' - [qExt(t,1:Lin) zeros(1,10)]'))))
%             logLikelihood = logLikelihood + logLikeliTemp;
%                         
%             logLikelihood = logLikelihood + logLikeliTemp;

    
    
    end
    BICMtx(thisItr,:) = temp4040;
    [~, indOfMinBIC] = min(temp4040);
    selcModel(thisItr) = indOfMinBIC;
    
    
%     regConst2 = 0.1;%0.00000001;
%     logLikelihoodEpi = 0;
%     logLikelihood = 0;
%     logLikelihoodSL = 0;
%     logLikelihood0 = 0;
%     timePointsUsed = 1;
%     for t = 1:numSamplingPoints-1
%         if(q(t,:) > 0 & q(t,:) < 1)
%             timePointsUsed = timePointsUsed + 1;
%             
%             
%             % making the V vetor for the Epi case
%             temp40 = zeros(numUnique2PtPairs, 1);
%             
%             count3 = 0;
%             for l2 = 1:Lin-1
%                 temp40(count3 + 1:count3 + Lin - l2) = temp40(count3 + 1:count3 + Lin - l2) + repmat(q(t,l2), Lin - l2, 1) + q(t, l2+1:end)';
%                 count3 = count3 + Lin - l2;
%             end
% 
%             vEpi = [ (1 - 2*q(t,:)');  (temp40-4*qExt(t,Lin+1:end)')];
%             
%             
% %             thisHmat = HmatAll(:,:,t);
% %             thisCmat = thisHmat(1:Lin,1:Lin);
% %             rankH = rank(thisHmat)
% %             rankC = rank(thisCmat)
% %             pause
%             
%     %         thisHmat
%     %         pause
%             
%             % all neutral SC
%             %den0 = det(thisCmat + regConst2*eye(Lin));
%             %vec0 = q(t+1,:)' - q(t,:)' - dT*(thisCmat + regConst2*eye(Lin))*zeros(5,1) - dT*muVal*(1 - 2*q(t,:)');
%             %logLikeliTemp0 = -Nin/(4*dT)*vec0'*(thisCmat + regConst2*eye(Lin))^-1*vec0 - 1/2*log(den0) + Lin/2*log(Nin/(4*pi*dT)) + (log(sum(abs(q(t+1,:)' - q(t,:)'))));
%             den0 = det(thisHmat + regConst2*eye(15));
%             vec0 = qExt(t+1,:)' - qExt(t,:)' - dT*thisHmat*zeros(15,1) - dT*muVal*vEpi;
%             logLikeliTemp0 = -Nin/(4*dT)*vec0'*(thisHmat + regConst2*eye(15))^-1*vec0 - 1/2*log(den0) + (Lin*(Lin + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(qExt(t+1,:)' - qExt(t,:)'))));
%             logLikelihood0 = logLikelihood0 + logLikeliTemp0;
% 
%             % SL SC
%             %thisCmatSL = diag(diag(thisCmat));
%             %denSL = det(thisCmatSL + regConst2*eye(Lin));
%             %vecSL = q(t+1,:)' - q(t,:)' - dT*(thisCmatSL + regConst2*eye(Lin))*sigmaEstOutUnLink - dT*muVal*(1 - 2*q(t,:)');
%             %logLikeliTempSL = -Nin/(4*dT)*vecSL'*(thisCmatSL + regConst2*eye(Lin))^-1*vecSL - 1/2*log(denSL) + Lin/2*log(Nin/(4*pi*dT)) + (log(sum(abs(q(t+1,:)' - q(t,:)'))));
%             thisHmatSL = diag([diag(thisCmat); zeros(10,1)]);
%             denSL = det(thisHmatSL + regConst2*eye(15));
%             vecSL = [qExt(t+1,1:Lin) zeros(1,10)]' - [qExt(t,1:Lin) zeros(1,10)]' - dT*thisHmatSL*[sigmaEstOutUnLink; zeros(10,1)] - dT*muVal*[vEpi(1:Lin); zeros(10,1)];
%             logLikeliTempSL = -Nin/(4*dT)*vecSL'*(thisHmatSL + regConst2*eye(15))^-1*vecSL - 1/2*log(denSL) + (Lin*(Lin + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs([qExt(t+1,1:Lin) zeros(1,10)]' - [qExt(t,1:Lin) zeros(1,10)]'))));
%             logLikelihoodSL = logLikelihoodSL + logLikeliTempSL;
%             
%             % MPL SC
%             den = det(thisCmat + regConst2*eye(Lin));
%             vec = q(t+1,:)' - q(t,:)' - dT*(thisCmat + regConst2*eye(Lin))*sigmaEstOutLink - dT*muVal*(1 - 2*q(t,:)');
%             logLikeliTemp = -Nin/(4*dT)*vec'*(thisCmat + regConst2*eye(Lin))^-1*vec - 1/2*log(den) + Lin/2*log(Nin/(4*pi*dT)) + (log(sum(abs(q(t+1,:)' - q(t,:)'))))
% %             den = det(thisHmat + regConst2*eye(15));
% %             vec = qExt(t+1,:)' - qExt(t,:)' - dT*thisHmat*[sigmaEstOutLink; zeros(10,1)] - dT*muVal*vEpi;
% %             logLikeliTemp = -Nin/(4*dT)*vec'*(thisHmat + regConst2*eye(15))^-1*vec - 1/2*log(den) + (Lin*(Lin + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(qExt(t+1,:)' - qExt(t,:)'))));
% %             thisHmatMPL = [thisCmat zeros(5,10); zeros(10,15)];
% %             den = det(thisHmatMPL + regConst2*eye(15));
% %             vec = [qExt(t+1,1:Lin) zeros(1,10)]' - [qExt(t,1:Lin) zeros(1,10)]' - dT*thisHmatMPL*[sigmaEstOutUnLink; zeros(10,1)] - dT*muVal*[vEpi(1:Lin); zeros(10,1)];
% %             logLikeliTemp = -Nin/(4*dT)*vec'*(thisHmatMPL + regConst2*eye(15))^-1*vec - 1/2*log(den) + (Lin*(Lin + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs([qExt(t+1,1:Lin) zeros(1,10)]' - [qExt(t,1:Lin) zeros(1,10)]'))))
%             logLikelihood = logLikelihood + logLikeliTemp;
%                         
%             logLikelihood = logLikelihood + logLikeliTemp;
%             
%             %MPLE SC
%             denEpi = det(thisHmat + regConst2*eye(15));
%             vecEpi = qExt(t+1,:)' - qExt(t,:)' - dT*thisHmat*sigmaEstOutLinkEpi - dT*muVal*vEpi;
%             logLikeliEpiTemp = -Nin/(4*dT)*vecEpi'*(thisHmat + regConst2*eye(15))^-1*vecEpi - 1/2*log(denEpi) + (Lin*(Lin + 1)/2)/2*log(Nin/(4*pi*dT)) + (log(sum(abs(qExt(t+1,:)' - qExt(t,:)'))))
%             logLikelihoodEpi = logLikelihoodEpi + logLikeliEpiTemp;
%             pause
%         end
%     end
%      
%     BICwithoutEpi0(thisItr) = (Lin*(Lin + 1)/2)*log(timePointsUsed*ng) - 2*logLikelihood0
%     BICwithoutEpiSL(thisItr) = (Lin*(Lin + 1)/2)*log(timePointsUsed*ng) - 2*logLikelihoodSL
%     BICwithoutEpi(thisItr) = (Lin*(Lin + 1)/2)*log(timePointsUsed*ng) - 2*logLikelihood % Lin*log(timePointsUsed*ng) - 2*logLikelihood;
%     BICwithEpi(thisItr) = (Lin*(Lin + 1)/2)*log(timePointsUsed*ng) - 2*logLikelihoodEpi
%     pause
%     [minValBIC(thisItr), indMinBIC(thisItr)] = min([BICwithoutEpi0(thisItr) BICwithoutEpiSL(thisItr) BICwithoutEpi(thisItr) BICwithEpi(thisItr)]);
%     %[minValBIC(thisItr), indMinBIC(thisItr)] = min([BICwithoutEpi0(thisItr)+1000000 BICwithoutEpiSL(thisItr) BICwithoutEpi(thisItr) BICwithEpi(thisItr)]);
    
    
    disp('done.')
    toc

    
    
    allEstEpi(thisItr,:) = sigmaEstOutLinkEpi';
    allEstEpiWithR(thisItr,:) = sigmaEstOutLinkEpiWithR';
    allEstNoMuEpi(thisItr,:) = sigmaEstOutLinkNoMuEpi';
    allEstEpiNoE(thisItr,:) = sigmaEstOutLinkEpiNoE';
    allEst(thisItr,:) = sigmaEstOutLink';
    
    
    
    % error in  s1, s1 and s12
    allNRMSEEpi(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi' - [diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2)./sum(abs([diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2));
    allNRMSE(thisItr) = sqrt(sum(abs([sigmaEstOutLink; zeros(length(perSiteSelctionAllEpiTerms),1)]' - [diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2)./sum(abs([diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2));
    
    allNRMSEEpiNoE(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpiNoE' - [diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2)./sum(abs([diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2));
    
    
    % save file
    if(loadDataOption == 1)        
        if(Tstart == 1)
            fileNameSave = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_reg1' regStr1 '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        else
             fileNameSave = [fileName(1:end-4) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_reg1' regStr1  '_reg2' regStr2 '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        end
    elseif(loadDataOption == 2)
        % in this case, only change the file extension
        
        temp60 = strfind(fileName, '_Set');
        
        fileNameSave = [fleName(1:temp60-1) '_reg1' regStr1 '_reg2' regStr2 fileName(temp60:end-4) '.mat'];
    end

    if(useStrongReg == 1)
        fileNameSave = [fileNameSave(1:end-4) '_useStrong.mat'];
    end
    timeWholeCodeMPL(thisItr) = toc;
    timeWholeCodeMPL(thisItr)

    if(saveFile == 1)
        disp('Saving file...')
        save([dirNameAnalysis fileNameSave])
    end
%     allNRMSEEpi(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2));
%     allNRMSE(thisItr) = sqrt(sum(abs([sigmaEstOutLink; 0]' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2));
%     
%     % error in only s1 and s2
%     allNRMSEEpi_Site1n2(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi(1:2)' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2));
%     allNRMSE_MPL_Site1n2(thisItr) = sqrt(sum(abs([sigmaEstOutLink]' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2));
%rref(Hmat)
%pause
%Hmat*dT% + Reg2Mat
%pause
end
