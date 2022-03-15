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
% code for 5-site simulation Figure 6

% Previous updated 10-Jan 2018
% Author: M Saqib Sohail
% 
clc

clear all
%close all

saveFile = 1;%1;

numStrainsInInitialPop = 20; % 5, 10, 30
%L = 50;
thisSet = 1064701;%1062001;%128001001;%158021001;%1561001;%10211;%5611002%7%5613%315%56%5611001%5601%5991%53900%6510;%6500%1990;%1991%1034;%7%58%7;
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

dTAll = [5 10 20 25 50];%%[5 10 20 30 50 100];%[5 10 15 20 30 50 75 100];
ngAll = [20 50 100 200];%[10 20 50 80 100 200];
for dd = 1:length(dTAll)
for nn = 1:length(ngAll)

dT = dTAll(dd)%10;
ng = ngAll(nn)%100;
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
end
end