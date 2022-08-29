
function [sumAijWithLink, q, v_EstOutLink, sumAijUnLink, sumEWithLink, ...
    sumFWithLink, qExt_tLast, qExt_tStart, vExt_EstOutLink, uExt_EstOutLinkWithR, ...
    dT, Reg1Mat, startTime, stopTime, Reg2Mat, loadDataOption, Tstart, ...
    fileName, dirNameAnalysis, Tend, ng, model, ...
    priorConst1, priorConst2, perSiteSelctionEpi] = run2siteMPLE(thisSet, thisItr)

loadDataOption = 1;
% 1 : load data from full trajectories information
% 2 : load data from sampled trajectories information (*.dat file)

priorConst1 = 1;%1/10;%1/N^2;
priorConst2 = 1;%1/10;%1/N^2;

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
elseif(thisSet >= 108701 && thisSet <= 118917)
    Tused = 150;
end
%Tend = Tused + Tstart - 1;
Tend = Tused + Tstart;    
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
    uExt_EstOutLinkWithR = 0;
    sigmaEstOutLinkEpiWithR = [0 0 0]';
end
