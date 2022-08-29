%% Run Maximum Path Likelihood (MPLin) method on GT data
% v2 is faster version. It calculates the int cov matrix in step 2 and does
% not save cov matrix at each time point

% Use different reg for epi terms
%%

% this code :
%      1. loads data (in .mat or .dat format)
%      2. calculates one and two point frequencies
%      3. Find estimates using MPL

% Last updated 10-Jan 2018
% Author: M Saqib Sohail

% 
clc
clear all
%close all

saveFile = 1;%1;

numStrainsInInitialPop = 20; % 5, 10, 30
%L = 50;
thisSet = 1064701;%1062001;%1064701;%1064001;%10211;%5611002%7%5613%315%56%5611001%5601%5991%53900%6510;%6500%1990;%1991%1034;%7%58%7;
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

timeWholeCodeIlling = zeros(1, numItr);
%%
numParam = Lin*(Lin+1)/2;
allEstEpi = zeros(numItr,numParam);
allEstNoMuEpi = zeros(numItr,numParam);
allEstEpiNoE = zeros(numItr,numParam);
allEst = zeros(numItr,Lin);
allEpiEstIlling = zeros(numItr,numParam);
allFval = zeros(numItr, 1);
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

        numGenos = 2^Lin;
        genotype = zeros(numGenos, Lin);
        colPatternRepTimes_temp = 2.^[0:(Lin-1)];
        ZeroOneRepTimes_temp = 2.^[(Lin-1):-1:0];
        for i = Lin:-1:1
            colPattern_temp = [zeros(ZeroOneRepTimes_temp(i), 1) ; 
                         ones(ZeroOneRepTimes_temp(i), 1)];                     
            colPattern = repmat(colPattern_temp, colPatternRepTimes_temp(i), 1);
            genotype(:,i) = colPattern;
        end

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

        N_k = -1*ones(numSamplingPoints, numGenos);
        n_i = -1*ones(numSamplingPoints, numGenos);
        qGeno = -1*ones(numSamplingPoints, numGenos);
        fprintf('Calculate 1 point probabilities...')
        for t = 1:numSamplingPoints

            thisSamplingTimePoint = samplingTimePoints(t);
            strainsThisTimePoint = masterTimeFreq{thisSamplingTimePoint}(:,1);
            freqStrainThisTimePoint = masterTimeFreq{thisSamplingTimePoint}(:,2);
            thisMSATemp = -1*ones(numSamplesPerGen, Lin);
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

            for l = 1:numGenos
                N_k(t,l) = numSamplesPerGenSelected;
                n_i(t,l) = sum(sum(thisMSA == genotype(l,:), 2) == Lin);
                qGeno(t, l) = n_i(t,l)/numSamplesPerGenSelected;
            end





        %    q(t,:) = n_i(t,:)./N_k(t); % normalize frequencies;
        end

        perSiteSelctionAllEpiTerms = [];
        for l = 1:Lin
            perSiteSelctionAllEpiTerms = [perSiteSelctionAllEpiTerms perSiteSelctionEpi(l,l+1:end)];
        end

        disp('done')

    end

    % plot tranjectory at each site
    if(plotFigs)
       figure
       plot(qGeno)
       xlabel('Generation')
       ylabel('Frequency')
       title(['Sampled genotype freq'])
    end

    % Define Hamming Distance and Mutation matrices for an L-locus
    % system 
    numGenoMaxPossible = 2^Lin;

    HD = -1*ones(numGenoMaxPossible, numGenoMaxPossible);
    for i = 1:numGenoMaxPossible
        geno1 = genotype(i,:);
        for j = 1:numGenoMaxPossible
            geno2 = genotype(j,:);
            HD(i,j) = sum(xor(geno1, geno2));
        end
    end

    M =  -1*ones(numGenoMaxPossible, numGenoMaxPossible);
    for i = 1:numGenoMaxPossible

        for j = 1:numGenoMaxPossible
            thisHD = HD(i,j);
            M(i,j) = muVal^thisHD*((1-muVal)^(Lin-thisHD));
        end
    end



    % 3. Compute likelihood    
    %--------------------------------------------------------------------------   
    
    selEpiVecGuess = [1:15]/100;
    %likelihood = eval_likelihood(selEpiVecGuess, Lin, numGenoMaxPossible, numSamplingPoints, dT, n_i, N_k, M, qGeno, genotype);

    numOptRuns = 1; % number of optimization runs
    selEpiVecGuessAll = zeros(Lin*(Lin+1)/2, numOptRuns);
    fvalAll = zeros(1, numOptRuns);
    SARuns = 5000;
    StallIterLimIn = 10000;
    TolFunIn = 0.0001/StallIterLimIn;
    dispStr = 'diagnose';%'diagnose';%'off';


    overAllLogLikelihood = @(selEpiVecGuessIn) eval_likelihood(selEpiVecGuessIn, Lin, numGenoMaxPossible, numSamplingPoints, dT, n_i, N_k, M, qGeno, genotype);

    lowerBound = -1*10^0*ones(1, Lin*(Lin+1)/2);
    upperBound = 1*10^0*ones(1, Lin*(Lin+1)/2);
    initialSigmaEstIn = 0.02*randn(1, Lin*(Lin+1)/2) - 0.01;
    optAlgo = 'SA';
    optionsSA = saoptimset('simulannealbnd');

    optionsSA = saoptimset(optionsSA,'MaxFunEval',SARuns,'TolFun',TolFunIn,'StallIterLim',StallIterLimIn,'Display',dispStr);
    [selEpiVecGuess, fval] = simulannealbnd(overAllLogLikelihood, initialSigmaEstIn, lowerBound, upperBound, optionsSA);
    allEpiEstIlling(thisItr,:) = selEpiVecGuess;
    allFval(thisItr) = fval;


    disp('done.')
 

    timeWholeCodeIlling(thisItr) = toc;
    timeWholeCodeIlling(thisItr)

    
    % save file
    if(loadDataOption == 1)        
        if(Tstart == 1)
            fileNameSave = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_Illing_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        else
             fileNameSave = [fileName(1:end-4) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_Illing_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
        end
    elseif(loadDataOption == 2)
        asd
    end

    if(saveFile == 1)
        disp('Saving file...')
        save([dirNameAnalysis fileNameSave])
    end

end
