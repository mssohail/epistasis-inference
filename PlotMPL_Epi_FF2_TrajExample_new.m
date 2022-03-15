%% Run Maximum Path Likelihood (MPLin) method on GT data
% v2 is faster version. It calculates the int cov matrix in step 2 and does
% not save cov matrix at each time point

% runs for new version of 2-site Epi test cases
% code that runs multiple sets at once


% this code :
%      1. loads data (in .mat or .dat format)
%      2. calculates one and two point frequencies
%      3. Find estimates using MPL

% Last updated 10-Jan 2018
% Author: M Saqib Sohail

% 
clc
clear all
close all

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

saveFigs = 1;

allSets = [201311 202211];
itrSets = [17 17]; %[5, 14, 15    6 3 15

figureDir = [pwd '\Figures\'];

for kk = 1:2
    %L = 50;
    thisSet = allSets(kk)%108711%108413;%108706;%10211;%5611002%7%5613%315%56%5611001%5601%5991%53900%6510;%6500%1990;%1991%1034;%7%58%7;
    numItr = 1000;%1000;%200%90%90%90;

    loadDataOption = 1;
    % 1 : load data from full trajectories information
    % 2 : load data from sampled trajectories information (*.dat file)

    priorConst1 = 1;
    priorConst2 = 1;%1/10;%1/N^2;

    runLinearInt = 0; % 1: run linear interpolation code, 0: skip it
    fileNameContainingDirPath = 'dirNames.txt';
    getSysParam;

    TgenStr = TgenStr(1:end);
    %TgenStr = num2str(str2double(TgenStr)*1000);
    dT = 10;%10;
    ng = 100;%100;
    Tstart = 1;
    Tused = 1000;%990;%300;%990;
    %Tend = Tused + Tstart - 1;
    Tend = Tused + Tstart;

    plotFigs = 0;

    timeWholeCodeMPL = zeros(1, numItr);


    allEstEpi = zeros(numItr,3);
    allEstEpiWithR = zeros(numItr,3);
    allEstNoMuEpi = zeros(numItr,3);
    allEstEpiNoE = zeros(numItr,3);
    allEst = zeros(numItr,2);

    for thisItr = itrSets(kk)
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
            % load data from full trajectories
            fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_selVal' num2str(selVal2) '_T' TgenStr '_' num2str(thisItr) '.mat'];

            load([dirNameData fileName], 'N', 'D', 'masterStrainList', 'masterTimeFreq', 'perSiteSelction', 'sitesUnderSelection', 'muVal', 'perSiteSelctionEpi_skeleton', 'perSiteSelctionEpi', 'currentNumOfCirStrainsInAllPop');


            masterStrainList = masterStrainList(:,1:Lin);
        elseif(loadDataOption == 2)
            % load data from sampled trajectories
            if(Tstart == 1)
                fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_selVal' num2str(delSelVal2) '_T' TgenStr '_' num2str(thisItr) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
            else
                fileName = ['WFsimEpi_Nmu' NmuInStr '_N' num2str(Nin) '_L' num2str(Lin) '_D' num2str(DAll) '_selVal' num2str(selVal1) '_selVal' num2str(delSelVal2) '_T' TgenStr '_' num2str(thisItr) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_linkDiff_Set' num2str(thisSet) '_itr' num2str(thisItr) '.dat'];
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
            q11Sum = zeros(numUnique2PtPairs, 1);
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

            saved4pt = -1*ones(1, numSamplingPoints);
            %randSelctIndAll = zeros(numSamplingPoints, numSamplesPerGenSelected);
            % this section calculates q_t for each site, for the synthetic protein, 
            % 2 allele per site case this is just the sum of occurances of ones
            fprintf('Calculate 1 and 2 point probabilites...')
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
                    tempq01 = (repmat((~thisMSALogical(:,l)),1, Lin).* thisMSALogical);
                    tempq10 = (repmat(thisMSALogical(:,l),1, Lin).* (~thisMSALogical));
                    tempq00 = (repmat((~thisMSALogical(:,l)),1, Lin).* (~thisMSALogical));

                    % sum for the lth row of the covariance matrix
                    %q11(l,:,t) = sum(tempq11)./numSamplesPerGenSelected;

                    q11Temp(l,:) = sum(tempq11)./numSamplesPerGenSelected;
                    q01Temp(l,:) = sum(tempq01)./numSamplesPerGenSelected;
                    q10Temp(l,:) = sum(tempq10)./numSamplesPerGenSelected;
                    q00Temp(l,:) = sum(tempq00)./numSamplesPerGenSelected;
                    q11(l,:,t) = q11Temp(l,:);
                    q01(l,:,t) = q01Temp(l,:);
                    q10(l,:,t) = q10Temp(l,:);
                    q00(l,:,t) = q00Temp(l,:);
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
                q11Sum = q11Sum + q11_jkTemp_partial'; %q11Sum + q11TempVec';  

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
            fprintf('Calculate 1 and 2 point probabilites...')
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


        sigmaEstOutLink = (sumAijWithLink*dT + priorConst1*eye(Lin))\(q(stopTime,:) - q(startTime,:) - v_EstOutLink')';
        sigmaEstOutLinkNoMu = (sumAijWithLink*dT + priorConst1*eye(Lin))\(q(stopTime,:) - q(startTime,:))';
        sigmaEstOutUnLink = (sumAijUnLink*dT + priorConst1*eye(Lin))\(q(stopTime,:) - q(startTime,:) - v_EstOutLink')';
        sigmaEstOutUnLinkNoMu = (sumAijUnLink*dT + priorConst1*eye(Lin))\(q(stopTime,:) - q(startTime,:))';

        Hmat = [sumAijWithLink sumEWithLink;
                sumEWithLink'  sumFWithLink];

        Reg2Mat = diag([priorConst1*ones(1,Lin) priorConst2*ones(1, numUnique2PtPairs)]);
        sigmaEstOutLinkEpi = (Hmat*dT + Reg2Mat)\(qExt_tLast - qExt_tStart - vExt_EstOutLink')';
        sigmaEstOutLinkNoMuEpi = inv(Hmat*dT + Reg2Mat)*(qExt_tLast - qExt_tStart)';

        % calc recomb term
        if(recombination == 1)
    %          recTwoPointsTemp = zeros(Lin*(Lin - 1)/2,1); % only epi trms will have this
    %          count1 = 1;
    %          for rr = 1:Lin-1
    %              mulVec = (1:Lin-rr)*(LinEff-1); % mulVec is magnitude of (i-j)
    %              numEntriesThis = length(mulVec);
    %              temprr2 = mulVec.*sumAijWithLink(rr, rr+1:Lin);
    %              recTwoPointsTemp(count1:count1+numEntriesThis-1) = temprr2';
    %              count1 = count1 + numEntriesThis;
    %          end
            temp33 = sumAijWithLink(1,2);%- q11Sum;% + recVal*recTwoPointsTemp + recVal*(LinEff-1)*q11Sum;
            uExt_EstOutLinkWithR = -dT*recVal*(LinEff-1)*[zeros(Lin, 1); temp33];

            sigmaEstOutLinkEpiWithR = (Hmat*dT + Reg2Mat)\(qExt_tLast - qExt_tStart - vExt_EstOutLink' - uExt_EstOutLinkWithR')';
        else
            sigmaEstOutLinkEpiWithR = [0 0 0]';
        end

    % [[perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]' sigmaEstOutLinkEpi sigmaEstOutLinkEpiWithR uExt_EstOutLinkWithR]
    %  pause   

        HmatNoE = [sumAijWithLink sumEWithLink*0;
                0*sumEWithLink'  sumFWithLink];
        sigmaEstOutLinkEpiNoE = (HmatNoE*dT + Reg2Mat)\(qExt_tLast - qExt_tStart - vExt_EstOutLink')';
        sigmaEstOutLinkNoMuEpiNoE = inv(HmatNoE*dT + Reg2Mat)*(qExt_tLast - qExt_tStart)';


        if(runLinearInt == 1)
            % LinearInter_Nsites_Dec2016;
            LinearInter_Nsites_Jan2017_1;

            startTime_LinUni = 1;
            stopTime_LinUni = size(q_LinUni,1);

            sumAijWithLink_LinUni = zeros(Lin,Lin);

            for t_LinUni = startTime_LinUni:stopTime_LinUni-1
                tempDiagEntries_LinUni = q_LinUni(t_LinUni,:).*(1 - q_LinUni(t_LinUni,:));
                temp200_LinUni = zeros(Lin,Lin);
                for thisSite = 1:Lin
                    % remember q11 is just the 2 point, i.e., q_ij, we still
                    % need ot subtract q_i*q_j from it
                    % this makes q_i*q_i q_i*q_j ... q_i*q_k ( will remove
                    % diagonal later)
                    temp100_LinUni = q_LinUni(t_LinUni,:).*q_LinUni(t_LinUni,thisSite);

                    % except for (i,i) entry, rest are q_ij - q_i*q_j
                    temp200_LinUni(thisSite,:) = q11_LinUni(thisSite, :, t_LinUni) - temp100_LinUni;
                end
                covMtxThisTime_LinUni = ~(eye(Lin)).*temp200_LinUni + diag(tempDiagEntries_LinUni);
                sumAijWithLink_LinUni = sumAijWithLink_LinUni + covMtxThisTime_LinUni;
            end
            tempTimeSumCovProcessig_LinUni = toc;

            v_EstOutLink_LinUni = muVal*(sum((ones(Lin,stopTime_LinUni-startTime_LinUni+1-1) - 2*q_LinUni(startTime_LinUni:stopTime_LinUni-1,:)'),2));
            sigmaEstOutLink_LinUni = (sumAijWithLink_LinUni*dT_LinUni + priorConst*eye(Lin))\(q_LinUni(stopTime_LinUni,:) - q_LinUni(startTime_LinUni,:) - v_EstOutLink_LinUni')';

        end    
        disp('done.')
        toc

        % save file
        if(loadDataOption == 1)        
            if(Tstart == 1)
                fileNameSave = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
            else
                 fileNameSave = [fileName(1:end-4) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
            end
        elseif(loadDataOption == 2)
            % in this case, only change the file extension
            fileNameSave = [fileName(1:end-4) '.mat'];
        end

        timeWholeCodeMPL(thisItr) = toc;
        timeWholeCodeMPL(thisItr)

        allEstEpi(thisItr,:) = sigmaEstOutLinkEpi';
        allEstEpiWithR(thisItr,:) = sigmaEstOutLinkEpiWithR';
        allEstNoMuEpi(thisItr,:) = sigmaEstOutLinkNoMuEpi';
        allEstEpiNoE(thisItr,:) = sigmaEstOutLinkEpiNoE';
        allEst(thisItr,:) = sigmaEstOutLink';    
    end

    twoPointFreqs{kk} = [squeeze(q00(1,2,:)) squeeze(q01(1,2,:)) squeeze(q10(1,2,:)) squeeze(q11(1,2,:))];
    covMtxCell{kk} = Hmat*dT + Reg2Mat;
    estimatesCell{kk} = sigmaEstOutLinkEpiWithR;
    perSiteSelctionEpiCell{kk} = perSiteSelctionEpi;
end


limits_data = [round(min(min(min(covMtxCell{1})), min(min(covMtxCell{2})))/10)*10 round(max(max(max(covMtxCell{1})), max(max(covMtxCell{2})))/10)*10];




%%
xDim = 13.4;
yDim = 8.6;
leftMargin = 1.25;
rightMargin = 0.25;
bottomMargin = 1;
bottomMarginp = 0.5;
topMargin = 0.5;

hgap1 = 0.25;%1.5;
hgap2 = 1.55;%2;
vgap = 1.5;

width1 = 4.5;%6.8;
width2 = 2.8;%4;
width3 = 2.8;%4;
height = 2.8;%4;

fig1 = figure('Units','centimeters', ...
                'Position', [10 3 xDim yDim]);%[10 3 14.1 13]);


ha(1) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin+vgap+height width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(2) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin+vgap+height width2 height], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(3) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1+hgap2+width2 bottomMargin+vgap+height width3 height], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(4) = axes('Units','centimeters', ...
                'Position',[leftMargin bottomMargin width1 height], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(5) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1 bottomMargin width2 height], ...
                'XTickLabel','', ...
                'YTickLabel','');

ha(6) = axes('Units','centimeters', ...
                'Position',[leftMargin+hgap1+width1+hgap2+width2 bottomMargin width3 height], ...
                'XTickLabel','', ...
                'YTickLabel','');
            
mapBlues = brewermap(8,'Blues');            
color_scheme1 = brewermap(100,'Blues');
color_scheme2 = brewermap(100,'Reds');
color_scheme3 = brewermap(100,'Greens');
color_scheme4 = brewermap(100,'Greys');
color_schemeHM1 = brewermap(64,'Greens');
white_scheme64 = repmat([ 1 1 1],64,1);
white_scheme = repmat([ 1 1 1],100,1);
alpha1 = 0.8;
alpha2 = 0.7;
color_scheme11 = (255*color_scheme1*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme21 = (255*color_scheme2*(alpha2) + 255*white_scheme*(1-alpha2))/255;
color_scheme31 = (255*color_scheme3*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_scheme41 = (255*color_scheme4*(alpha1) + 255*white_scheme*(1-alpha1))/255;
color_schemeHM = (255*color_schemeHM1*(alpha1) + 255*white_scheme64*(1-alpha1))/255;
color_bar_visibility = 'off';


for kk = 1:2
    twoPoints = twoPointFreqs{kk};
    
    axes(ha((kk-1)*3+1))
    plot(0:10:1000, twoPoints(:,1),'-', 'color', color_scheme41(70,:), 'LineWidth', 1)
    hold on
    plot(0:10:1000, twoPoints(:,2),'-', 'color', color_scheme11(90,:), 'LineWidth', 1)
    plot(0:10:1000, twoPoints(:,3),'-', 'color', color_scheme21(70,:), 'LineWidth', 1)
    plot(0:10:1000, twoPoints(:,4),'-', 'color', color_scheme31(70,:), 'LineWidth', 1)
    ylabel('Frequency')
    xlabel('Generation')
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
      'XTick', [0 250 500 750 1000], ...%       
      'YTick', [0 0.25 0.5 0.75 1], ...'XTickLabel', {'s_1', 's_2', 's_{12}'}, ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
      'LineWidth', 0.5)
    if(kk == 1)
        leg = legend('00', '01', '10', '11');
        set(leg,'color','none');
        set(leg, 'Edgecolor','none');
        set(leg, 'units', 'centimeter', 'Position', [3.8 6.15 1 1])
        leg.Title.String = 'Genotypes';
    end
    
    
    axes(ha((kk-1)*3+2))
    
    
    hmap1 = heatmap(1:3, 1:3, covMtxCell{kk}, 'ColorBarVisible', color_bar_visibility);
    emptyCell = cell(3,1);
    emptyCell(1) = {' '};
    emptyCell(2) = {' '};
    emptyCell(3) = {' '};
    hmap1.YDisplayLabels = emptyCell;
    hmap1.XDisplayLabels = emptyCell;
    hmap1.FontSize = 8;
    hmap1.FontName = 'Arial';
    hmap1.Colormap = color_schemeHM;
    hmap1 = set(gca,'ColorLimits',limits_data);
    
    
    perSiteSelctionEpiTemp = perSiteSelctionEpiCell{kk};
    trueSel(1) = perSiteSelctionEpiTemp(1,1);
    trueSel(2) = perSiteSelctionEpiTemp(2,2);
    trueSel(3) = perSiteSelctionEpiTemp(1,2);
    
    
    invMtx = inv(covMtxCell{kk});
    stdErr = 1.96*sqrt(diag(invMtx)/1000);
    estimates = estimatesCell{kk};
    
    axes(ha((kk-1)*3+3))
    plot(0.7:0.1:1.3, trueSel(1)*ones(1, 7), 'color', color_scheme41(50,:), 'LineWidth', 2)
    hold on
    plot(1.7:0.1:2.3, trueSel(2)*ones(1, 7), 'color', color_scheme41(50,:), 'LineWidth', 2)
    plot(2.7:0.1:3.3, trueSel(3)*ones(1, 7), 'color', color_scheme41(50,:), 'LineWidth', 2)
    
    step1 = stdErr(1)/10;
    temp11 = (estimates(1)-stdErr(1)):step1:(estimates(1)+stdErr(1));
    temp12 = length(temp11);
    step2 = stdErr(2)/10;
    temp21 = (estimates(2)-stdErr(2)):step2:(estimates(2)+stdErr(2));
    temp22 = length(temp21);
    step3 = stdErr(3)/10;
    temp31 = (estimates(3)-stdErr(3)):step3:(estimates(3)+stdErr(3));
    temp32 = length(temp31);
    plot(ones(1, temp12), temp11, 'LineWidth', 1, 'color', color_scheme41(95,:))
    hold on
    plot(2*ones(1, temp22), temp21, 'LineWidth', 1, 'color', color_scheme41(95,:))
    plot(3*ones(1, temp32), temp31, 'LineWidth', 1, 'color', color_scheme41(95,:))
    plot(1:1:3, estimates, 'o', 'MarkerFaceColor', color_scheme21(60,:), 'MarkerEdgeColor', color_scheme41(95,:), 'MarkerSize', 4)
    
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
      'XTick', [1:1:3], ...%       
      'YTick', [-0.05:0.025:0.05], ...
      'XTickLabel', {'s_1', 's_2', 's_{12}'}, ... 'YLim', [0.5 1.001], ...'XLim', [0.5 sum(indSelected_logical) + 0.5], ...'YTick'       , 0:0.1:1, ...'XTick'       , 1.5:0.5:2.5, ...
      'LineWidth', 0.5)
    ylabel('Estimate')
    axis([0.5 3.5 -0.05 0.05])
    
    
    
end

dimDummy = [0.1 0.1 0.1 0.1]; % dummy position (in normalized units)

tb = annotation('textbox',dimDummy,'String','A','FitBoxToText','on')
tb.Units = 'centimeter';
tb.Position = [-0.1 8.2 0.7 0.5];
tb.LineStyle = 'none';
tb.FontWeight = 'bold';
tb.FontName = 'Arial';
tb.FontSize = 12;

tc = annotation('textbox',dimDummy,'String','B','FitBoxToText','on')
tc.Units = 'centimeter';
tc.Position = [-0.1 3.9 0.7 0.5];
tc.LineStyle = 'none';
tc.FontWeight = 'bold';
tc.FontName = 'Arial';
tc.FontSize = 12;

ta1 = annotation('textbox',dimDummy,'String','Integrated','FitBoxToText','on')
ta1.Units = 'centimeter';
ta1.Position = [6.5 3.2 3 02];
ta1.LineStyle = 'none';
ta1.FontName = 'Arial';
ta1.FontSize = 8;

ta2 = annotation('textbox',dimDummy,'String','covariance matrix','FitBoxToText','on')
ta2.Units = 'centimeter';
ta2.Position = [6 2.9 3 02];
ta2.LineStyle = 'none';
ta2.FontName = 'Arial';
ta2.FontSize = 8;

tb1 = annotation('textbox',dimDummy,'String','Integrated','FitBoxToText','on')
tb1.Units = 'centimeter';
tb1.Position = [6.5 -1 3 02];
tb1.LineStyle = 'none';
tb1.FontName = 'Arial';
tb1.FontSize = 8;

tb2 = annotation('textbox',dimDummy,'String','covariance matrix','FitBoxToText','on')
tb2.Units = 'centimeter';
tb2.Position = [6 -1.3 3 02];
tb2.LineStyle = 'none';
tb2.FontName = 'Arial';
tb2.FontSize = 8;
if(exist(figureDir, 'dir') == 0)
    mkdir(figureDir)        
end

if(saveFigs == 1)
    figname = [figureDir 'Figure2_Example_trajs_2'];
    if(exist([figname '.jpg'], 'file') == 2)
        delete(figname)        
    end
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 xDim yDim])%,[0 0 8 6.45])% ,[0 0 8 6])
    %set(gcf, 'renderer', 'painters');
    print(figname, '-djpeg','-r400')
%     pause(0.5)
%     print(figname, '-dpng','-r400')
end
