%% Run Maximum Path Likelihood (MPLin) method on GT data
% v2 is faster version. It calculates the int cov matrix in step 2 and does
% not save cov matrix at each time point

% runs for new version of 2-site Epi test cases
% code that runs multiple sets at once
%%

% this code :
%      1. loads data (in .mat or .dat format)
%      2. calculates one and two point frequencies
%      3. Find estimates using MPL

% Last updated 03-Jan-2021
% code for 2-site simulations. Analysis for Figure 1 (8 cases), Figure 2, 
% Figure 3 and S1 (data variability)

% previous update 10-Jan 2018
% Author: M Saqib Sohail



% 
clc

clear all
close all

saveFile = 1;%1;
analyzeDataForFigure = 1; % 3;#

if(analyzeDataForFigure == 1)
    allSets = [201111 201211 201311 201411 202111 202211 202311 202411];% 102111 102211 102311 102411]%118711:118717%[108701:108707 108711:108717]
elseif(analyzeDataForFigure == 3)
    allSets = [118111:118117 118311:118317 118411:118417 118511:118517 118711:118717];
else
    disp('Error: Not a vald choice for variable analyzeDataForFigure')
end
for kk = 1:length(allSets)

    thisSet = allSets(kk)%108711%108413;%108706;%10211;%5611002%7%5613%315%56%5611001%5601%5991%53900%6510;%6500%1990;%1991%1034;%7%58%7;
    numItr = 1000;%1000;%200%90%90%90;

    loadDataOption = 1;
    % 1 : load data from full trajectories information
    % 2 : load data from sampled trajectories information (*.dat file)

    priorConst1 = 1;%1/10;%1/N^2;
    priorConst2 = 2;%1/10;%1/N^2;

    runLinearInt = 0; % 1: run linear interpolation code, 0: skip it
    fileNameContainingDirPath = 'dirNames.txt';
    getSysParam;

    TgenStr = TgenStr(1:end);
    %TgenStr = num2str(str2double(TgenStr)*1000);
    dT = 10;%10;
    ng = 100;%100;
    Tstart = 1;
    if(thisSet >= 201111 && thisSet <= 208411)
        Tused = 1000;%990;%300;%990;
    elseif(thisSet >= 108701 && thisSet <= 118717)
        Tused = 150;
    end
    %Tend = Tused + Tstart - 1;
    Tend = Tused + Tstart;

    plotFigs = 0;

    timeWholeCodeMPL = zeros(1, numItr);
    %%

    allEstEpi = zeros(numItr,3);
    allEstEpiWithR = zeros(numItr,3);
    allEstNoMuEpi = zeros(numItr,3);
    allEstEpiNoE = zeros(numItr,3);
    allEst = zeros(numItr,2);
    allNRMSEEpi = zeros(1, numItr);
    allNRMSEEpiWithR = zeros(1, numItr);
    allNRMSEEpiNoE = zeros(1, numItr);
    allNRMSEEpi_Site1n2 = zeros(1, numItr);
    allNRMSE = zeros(1, numItr);
    allNRMSE_MPL_Site1n2 = zeros(1, numItr);

    allNRMSEEpi_s1 = zeros(1, numItr);
    allNRMSEEpi_s2 = zeros(1, numItr);
    allNRMSEEpi_s12 = zeros(1, numItr);
    allNRMSEEpi_s1n2 = zeros(1, numItr);
    allNRMSEEpi_s1n12 = zeros(1, numItr);
    allNRMSEEpi_s2n12 = zeros(1, numItr);
    allNRMSEEpi_s1n2n12 = zeros(1, numItr);

    for thisItr = 1:numItr
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

        Reg1Mat = priorConst1*eye(Lin);
        Reg2Mat = diag([priorConst1*ones(1,Lin) priorConst2*ones(1, numUnique2PtPairs)]);

        sigmaEstOutLink = (sumAijWithLink*dT + Reg1Mat)\(q(stopTime,:) - q(startTime,:) - v_EstOutLink')';
        sigmaEstOutLinkNoMu = (sumAijWithLink*dT + Reg1Mat)\(q(stopTime,:) - q(startTime,:))';
        sigmaEstOutUnLink = (sumAijUnLink*dT + Reg1Mat)\(q(stopTime,:) - q(startTime,:) - v_EstOutLink')';
        sigmaEstOutUnLinkNoMu = (sumAijUnLink*dT + Reg1Mat)\(q(stopTime,:) - q(startTime,:))';

        Hmat = [sumAijWithLink sumEWithLink;
                sumEWithLink'  sumFWithLink];
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

    %     if(saveFile == 1)
    %         disp('Saving file...')
    %         save([dirNameAnalysis fileNameSave])
    %     end

        allEstEpi(thisItr,:) = sigmaEstOutLinkEpi';
        allEstEpiWithR(thisItr,:) = sigmaEstOutLinkEpiWithR';
        allEstNoMuEpi(thisItr,:) = sigmaEstOutLinkNoMuEpi';
        allEstEpiNoE(thisItr,:) = sigmaEstOutLinkEpiNoE';
        allEst(thisItr,:) = sigmaEstOutLink';



        % error in  s1, s1 and s12
        allNRMSEEpi(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi' - [diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2)./sum(abs([diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2));
        allNRMSEEpiWithR(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpiWithR' - [diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2)./sum(abs([diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2));
        allNRMSE(thisItr) = sqrt(sum(abs([sigmaEstOutLink; zeros(length(perSiteSelctionAllEpiTerms),1)]' - [diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2)./sum(abs([diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2));
        allNRMSEEpiNoE(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpiNoE' - [diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2)./sum(abs([diag(perSiteSelctionEpi)' perSiteSelctionAllEpiTerms]).^2));

        allNRMSEEpi_s1(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi(1) - perSiteSelctionEpi(1,1)).^2)./sum(abs(perSiteSelctionEpi(1,1)).^2));
        allNRMSEEpi_s2(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi(2) - perSiteSelctionEpi(2,2)).^2)./sum(abs(perSiteSelctionEpi(2,2)).^2));
        allNRMSEEpi_s12(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi(3) - perSiteSelctionEpi(1,2)).^2)./sum(abs(perSiteSelctionEpi(1,2)).^2));
        allNRMSEEpi_s1n2(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi(1:2)' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2));
        allNRMSEEpi_s1n12(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi([1 3])' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(1,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(1,2)]).^2));
        allNRMSEEpi_s2n12(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi([2 3])' - [perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2)./sum(abs([perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2));
        allNRMSEEpi_s1n2n12(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi([1 2 3])' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2));


    %     allNRMSEEpi(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2));
    %     allNRMSE(thisItr) = sqrt(sum(abs([sigmaEstOutLink; 0]' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2) perSiteSelctionEpi(1,2)]).^2));
    %     
    %     % error in only s1 and s2
    %     allNRMSEEpi_Site1n2(thisItr) = sqrt(sum(abs(sigmaEstOutLinkEpi(1:2)' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2));
    %     allNRMSE_MPL_Site1n2(thisItr) = sqrt(sum(abs([sigmaEstOutLink]' - [perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2)./sum(abs([perSiteSelctionEpi(1,1) perSiteSelctionEpi(2,2)]).^2));
    end

    if(recombination == 1)
        PlotMPL_Epi_2sites_recomb_funcBoxplot;
    else
        PlotMPL_Epi_2sites_funcBoxplot;
    end

    if(Tstart == 1)
        fileNameSaveCombined = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_reg1' num2str(priorConst1)  '_reg2' num2str(priorConst2) '_Set' num2str(thisSet) '_itr1_' num2str(thisItr) fileName(end-3:end)];
    else
         fileNameSaveCombined = [fileName(1:end-4) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_reg1' num2str(priorConst1)  '_reg2' num2str(priorConst2) '_Set' num2str(thisSet) '_itr1_' num2str(thisItr) fileName(end-3:end)];
    end

    if(saveFile == 1)
        disp('Saving file...')
        save([dirNameAnalysis fileNameSaveCombined])
    end


end

color_scheme_npg = [...
    0.9020    0.2941    0.2078; ...
    0.3020    0.7333    0.8353; ...
         0    0.6275    0.5294; ...
    0.2353    0.3294    0.5333; ...
    0.9529    0.6078    0.4980; ...
    0.5176    0.5686    0.7059; ...
    0.5686    0.8196    0.7608; ...
    0.8627         0         0; ...
    0.4941    0.3804    0.2824; ...
    0.6902    0.6118    0.5216 ];
color_scheme_npg  = [color_scheme_npg ; color_scheme_npg ];

chosenColors = color_scheme_npg([4 5 6],:);
figure
for l = 1:2
    plot(Tstart:dT:Tend, q(:,l), 'color', chosenColors(l,:), 'lineWidth', 1.5)
    hold on
end
plot(Tstart:dT:Tend, squeeze(q11(1,2,:)), 'o', 'color', chosenColors(3,:), 'lineWidth', 1.5)
xlabel('Generations')
ylabel('Freqeuncy')
axis([1 160 0 1])

leg = legend('q_1', 'q_2', 'q_{12}', 'location', 'SouthEast');
%leg = legend('MPLE', 'MPL', 'SL', 'location', 'NorthEast');

set(leg,'color','none');
set(leg, 'Edgecolor','none');

saveFigs = 0
% save Matlab Figure
if(saveFigs)
   %figname = [dirNameScriptFile chosenSlash 'Fig_Traj_plot6cases_' num2str(floor(thisSet/10)) '_GTWith' baseCase'];
   figname = ['Fig_Traj_plot6cases_' num2str(thisSet) '_GTWith' baseCase];
   set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 13])%,[0 0 8 6.45])% ,[0 0 8 6])
   set(gcf, 'renderer', 'painters');
   print(figname, '-dpng','-r400')
end