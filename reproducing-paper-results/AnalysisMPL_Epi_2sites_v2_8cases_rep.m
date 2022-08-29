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
%close all

saveFile = 1;%1;
analyzeDataForFigure = 3%1; % 3;#



%[118716 118516 118316 118916 118216 118016;
% 1: only 00 01
% 2: only 00 10
% 3: only 00 11
% 4: only 00 01 10
% 5: only 00 01 11
% 6: only 00 10 11
% 7: only 00 01 10 11             
    
bothSets  = [118713 118714];
plotFigs = 0;


numItr = 1000;%1000;%200%90%90%90;
timeWholeCodeMPL = zeros(1, numItr);
%%

allEstEpi_rep = zeros(numItr,3);
allEstEpiWithR_rep = zeros(numItr,3);
allEstNoMuEpi_rep = zeros(numItr,3);
allEstEpiNoE_rep = zeros(numItr,3);
allEst_rep = zeros(numItr,2);

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

    thisSet_1 = bothSets(1);
    thisSet_2 = bothSets(2);
    


    
    [sumAijWithLink_1, q_1, v_EstOutLink_1, sumAijUnLink_1, sumEWithLink_1, ...
        sumFWithLink_1, qExt_tLast_1, qExt_tStart_1, vExt_EstOutLink_1, uExt_EstOutLinkWithR_1] = run2siteMPLE(thisSet_1, thisItr);
    [sumAijWithLink_2, q_2, v_EstOutLink_2, sumAijUnLink_2, sumEWithLink_2, ...
        sumFWithLink_2, qExt_tLast_2, qExt_tStart_2, vExt_EstOutLink_2, uExt_EstOutLinkWithR_2, ...
        dT, Reg1Mat, startTime, stopTime, Reg2Mat, loadDataOption, Tstart, ...
        fileName, dirNameAnalysis,  Tend, ng, model, ...
        priorConst1, priorConst2, perSiteSelctionEpi] = run2siteMPLE(thisSet_2, thisItr);

    sigmaEstOutLink_rep = (sumAijWithLink_1*dT + sumAijWithLink_2*dT + Reg1Mat)\((q_1(stopTime,:) - q_1(startTime,:) - v_EstOutLink_1') + (q_2(stopTime,:) - q_2(startTime,:) - v_EstOutLink_2'))';
    sigmaEstOutLinkNoMu_rep = (sumAijWithLink_1*dT + sumAijWithLink_2*dT + Reg1Mat)\((q_1(stopTime,:) - q_1(startTime,:)) + (q_2(stopTime,:) - q_2(startTime,:)))';
    sigmaEstOutUnLink_rep = (sumAijUnLink_1*dT + sumAijUnLink_2*dT + Reg1Mat)\((q_1(stopTime,:) - q_1(startTime,:) - v_EstOutLink_1') + (q_2(stopTime,:) - q_2(startTime,:) - v_EstOutLink_2'))';
    sigmaEstOutUnLinkNoMu_rep = (sumAijUnLink_1*dT + sumAijUnLink_2*dT + Reg1Mat)\(q_1(stopTime,:) - q_2(startTime,:))';

    Hmat_1 = [sumAijWithLink_1 sumEWithLink_1;
              sumEWithLink_1'  sumFWithLink_1];
    Hmat_2 = [sumAijWithLink_2 sumEWithLink_2;
              sumEWithLink_2'  sumFWithLink_2];
    sigmaEstOutLinkEpi_rep = (Hmat_1*dT + Hmat_2*dT + Reg2Mat)\((qExt_tLast_1 - qExt_tStart_1 - vExt_EstOutLink_1') + (qExt_tLast_2 - qExt_tStart_2 - vExt_EstOutLink_2'))';
    sigmaEstOutLinkNoMuEpi_rep = inv(Hmat_1*dT + Hmat_2*dT + Reg2Mat)*((qExt_tLast_1 - qExt_tStart_1) + (qExt_tLast_2 - qExt_tStart_2))';

    
    
    



    
    disp('done.')
    toc

%     % save file
%     if(loadDataOption == 1)        
%         if(Tstart == 1)
%             fileNameSave = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
%         else
%              fileNameSave = [fileName(1:end-4) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_Set' num2str(thisSet) '_itr' num2str(thisItr) fileName(end-3:end)];
%         end
%     elseif(loadDataOption == 2)
%         % in this case, only change the file extension
%         fileNameSave = [fileName(1:end-4) '.mat'];
%     end

    timeWholeCodeMPL(thisItr) = toc;
    timeWholeCodeMPL(thisItr)

%     if(saveFile == 1)
%         disp('Saving file...')
%         save([dirNameAnalysis fileNameSave])
%     end

    allEstEpi_rep(thisItr,:) = sigmaEstOutLinkEpi_rep';
    %allEstEpiWithR_rep(thisItr,:) = sigmaEstOutLinkEpiWithR_rep';
    allEstNoMuEpi_rep(thisItr,:) = sigmaEstOutLinkNoMuEpi_rep';
    allEst_rep(thisItr,:) = sigmaEstOutLink_rep';


end

%%


if(Tstart == 1)
    fileNameSaveCombined = [fileName(1:end-4) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_reg1' num2str(priorConst1)  '_reg2' num2str(priorConst2) '_Set' num2str(thisSet_1) '_Set' num2str(thisSet_2) '_itr1_' num2str(thisItr) '_rep' fileName(end-3:end)];
else
     fileNameSaveCombined = [fileName(1:end-4) '_Tstart' num2str(Tstart) '_Tend' num2str(Tend) '_dts' num2str(dT) '_ng' num2str(ng) '_' model '_reg1' num2str(priorConst1)  '_reg2' num2str(priorConst2) '_Set' num2str(thisSet_1) '_Set' num2str(thisSet_2) '_itr1_' num2str(thisItr) '_rep' fileName(end-3:end)];
end

if(saveFile == 1)
    disp('Saving file...')
    save([dirNameAnalysis fileNameSaveCombined])
end
    
    
    
    
    
    
    
    
    
    
  