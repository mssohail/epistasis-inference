function analysisMPL_Epi_v2_outerFunction(priorConstSC, priorConstEpi, ...
    recombProb, refernceSequence, dataDirNameMain, analysisDirNameMain, ...
    resultsDirNameMain, fileNameContainingDirPathCell, fileNameContainingMutProb, lengthAlignedSeq)


% specify the number of NT that 'can' occur in the provided fasta files.
% Default value is 5 (ACGT-) 
numNT = 5; 



str1 = fileNameContainingDirPathCell{1};
% preprocessing Step 0
preprocessingStep_0_outerFunction(dataDirNameMain, analysisDirNameMain, resultsDirNameMain, str1)
% preprocessing Step 1
preprocessingSte_1_outerFunction(fileNameContainingDirPathCell, fileNameContainingMutProb, numNT, lengthAlignedSeq)

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

FLAG_MarkAccessibility = true; % KEEP this FALSE for the time being, Accessibility code needs to be checked

FLAG_SaveIntCovMtx = false; % SET: will save Integrated Covariance matrix (for debugging only)
FLAG_useFreqEntry = true;
FLAG_troubleShoot = false; % SET: saves SelEstNoMu and SelEstSLNoMu
FLAG_Epi = true; % SET: use MPL with epistasis, UNSET: MPL with epistasis not used
FLAG_UserProvidedRefSeq = true; % SET: user provides reference sequence in ACGT- form
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


fileNamesListThisDir = findFileNamesWithGivenText(dirNameStr1Files, fileNameContainingDirPathCell);
numPat = length(fileNamesListThisDir);
% -------------------------------------------------------------------------
if(numPat == 0)
    disp('NumPat = 0. Check initialization settings and run again.')
end

% ========================== BEGIN PROCESSING =============================


fileNameContainingDirPath = [dirNameStr1Files fileNamesListThisDir{1}];

FLAG_Skip = false;
if(FLAG_Skip == false)

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
