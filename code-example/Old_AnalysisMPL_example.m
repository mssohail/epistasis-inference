clc
clear all
close all
warning off

    %========================== INITIALIZATION ================================

    % -------------------------- User specified -------------------------------
%     thisSet = allSets(asets);
%     dTStep = dTAll(asets);
%     ng = ngAll(asets);
%     Tused = TusedAll(asets);
%     numStrainsInInitialPop = allStrains(asets);%20; % number of strains in the initial population


    % chose convention 1: Ito, 2: Stratonovich, 3: Linear interpolation
    % Stratonovich has increased robustness to sampling effects than Ito
    % Linear interpolation has most increased robustness to sampling effects
    setConvention = 1;

    priorConstSC = 1; % this is the strength of the SC regularization term
    priorConstEpi = 1; % this is the strength of the Epi regularization term

    thisGenomicSegStartInd = 1; % this is the starting location of the first NT of the protein in the whole genome
    thisGenomicSegStopInd = 5; % this is the ending location of the last NT of the protein in the whole genome

    % this file contains the NT-to-NT mutation probability. It must be located
    % in the folder .../MPL Pipeline/Data_Misc/MutationProbabilities/
    fileNameContainingMutProb = 'MutProb_SyntheticData_1eminus4.txt';
    recombProb = 1e-4; % recombination probability

    FLAG_UserProvidedRefSeq = true; % SET: user provides reference sequence in ACGT- form
    refernceSequence = repmat('A', 1, 5);
    FLAG_binaryApprox = true; % SET: use binary approximation (only binary approximation wroks currently)
    numNT = 5; % specify the number of NT that 'can' occur in the provided fasta files. Default value is 5 (ACGT-)

    FLAG_MarkAccessibility = false; % KEEP this FALSE for the time being, Accessibility code needs to be checked

    FLAG_SaveIntCovMtx = false; % SET: will save Integrated Covariance matrix (for debugging only)
    FLAG_useFreqEntry = true;
    FLAG_troubleShoot = false; % SET: saves SelEstNoMu and SelEstSLNoMu
    FLAG_Epi = true; % SET: use MPL with epistasis, UNSET: MPL with epistasis not used

    textCell{1} = 'dirNamesExample';%['dirNamesSet' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dTStep) '_Tused' num2str(Tused) '_initStr' num2str(numStrainsInInitialPop) '_' ];
    % -------------------------------------------------------------------------


    % ------------------------- AUTO INITIALIZATION ---------------------------
    % NO USER INPUT REQUIRED
    if(setConvention == 1)
        FLAG_stratonovich = false;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
        FLAG_linearInt = false;
    elseif(setConvention == 2)
        FLAG_stratonovich = true;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
        FLAG_linearInt = false;
    elseif(setConvention == 3)
        FLAG_stratonovich = false;%true;%true; % SET: stratonovich convention, UNSET: Ito convention
        FLAG_linearInt = true;
    end

    % this file will contain the names of .fasta files to analyze suing the
    % AnalysisMPL_shortRead code. This file will be generated autotomatically
    % in preprocessingStep1. Here we just need to specify the name. 
    fileNameFastaFilesWithHU = 'fastaFilesHU.txt'; 
    meFastaFileToAnalyze = 'fastaFilesToAnalyze.txt'; % right now, these fasta files need to be generated on laptop

    FLAG_firstSeqIsRef = true; % set: 1st sequence of every fasta file is reference sequence
    mainDir = pwd;
    if(ispc)
        chosenSlash = '\';
    elseif(isunix)
        chosenSlash = '/';
    else
        disp('Error: system is not unix and not PC...')
        pause
    end
    dirNameTemp123 = 'dirNameFiles';
    dirNameStr1Files = [mainDir chosenSlash 'Data_Misc' chosenSlash dirNameTemp123 chosenSlash];


    fileNamesListThisDir = findFileNamesWithGivenText(dirNameStr1Files, textCell);
    numPat = length(fileNamesListThisDir);
    % -------------------------------------------------------------------------
    if(numPat == 0)
        disp('NumPat = 0. Check initialization settings and run again.')
    end

    % ========================== BEGIN PROCESSING =============================

%     for pat = 1:numPat

        fileNameContainingDirPath = [dirNameStr1Files fileNamesListThisDir{1}];
%         indOfDash = strfind(fileNameContainingDirPath, '_');
%         indOfDot = strfind(fileNameContainingDirPath, '.');
%         patID = fileNameContainingDirPath(indOfDash(end-1)+1:indOfDash(end)-1);
%         thisProt = fileNameContainingDirPath(indOfDash(end)+1:indOfDot(end)-1);
% 
%         disp('-----------------------------------------------------------------')
%         disp(' ')
%         disp(['Patient: ' patID])
%         disp(['Protein: ' thisProt])
        FLAG_Skip = false;
%         if(strcmp(patID, 'p3') && strcmp(thisProt, 'p6'))
%             FLAG_Skip = true;
%         end
        if(FLAG_Skip == false)
            %analysisStep1_v2(fileNameContainingDirPath, priorConstSC, FLAG_stratonovich, FLAG_MarkAccessibility, FLAG_UserProvidedRefSeq, FLAG_SaveIntCovMtx, FLAG_useFreqEntry, FLAG_troubleShoot, FLAG_linearInt);
            priorConst = [priorConstSC;
                          priorConstEpi;
                          recombProb];
            FLAG_vector = [FLAG_stratonovich;
                           FLAG_MarkAccessibility;
                           FLAG_UserProvidedRefSeq;
                           FLAG_SaveIntCovMtx;
                           FLAG_useFreqEntry;
                           FLAG_troubleShoot;
                           FLAG_linearInt;
                           FLAG_Epi];

            analysisStep1_v2(fileNameContainingDirPath, priorConst, FLAG_vector, refernceSequence);
        end
%     end
% end
