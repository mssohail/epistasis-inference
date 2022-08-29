function analysisMPL_Epi_v2_rep_outerFunction(priorConstSC, priorConstEpi, ...
    recombProb, refernceSequence, dataDirNameMainCell, analysisDirNameMainCell, ...
    resultsDirNameMainCell, fileNameContainingDirPathCell, fileNameContainingMutProb, ...
    lengthAlignedSeq, numRep, dirNameAnalysisRep)


% specify the number of NT that 'can' occur in the provided fasta files.
% Default value is 5 (ACGT-) 
numNT = 5; 


for i = 1:numRep

    str1 = fileNameContainingDirPathCell{i};
    % preprocessing Step 0
    preprocessingStep_0_outerFunction(dataDirNameMainCell{i}, analysisDirNameMainCell{i}, resultsDirNameMainCell{i}, str1)
    % preprocessing Step 1
    preprocessingSte_1_outerFunction({fileNameContainingDirPathCell{i}}, fileNameContainingMutProb, numNT, lengthAlignedSeq)

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
    FLAG_RepComb = true;
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


    fileNamesListThisDir = findFileNamesWithGivenText(dirNameStr1Files, {fileNameContainingDirPathCell{i}});
    numPat = length(fileNamesListThisDir);
    % -------------------------------------------------------------------------
    if(numPat == 0)
        disp('Variable fileNameContainingDirPathCell is empty. Check initialization settings and run again.')
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
                       FLAG_Epi; 
                       FLAG_RepComb];

        analysisStep1_v2(fileNameContainingDirPath, priorConst, FLAG_vector, refernceSequence);
    end
end




% estimate from replicate combining

if(exist(dirNameAnalysisRep, 'dir') == 0)
    mkdir(dirNameAnalysisRep)        
end


intCovMtxEpiAll = 0;
numerEpiAll = 0;
for i = 1:numRep
    load([analysisDirNameMainCell{i} 'Estimates' chosenSlash 'tempDataFileforRepComb.mat'], ...
        'intCovMtxEpi', 'numerEpi', 'regMtxEpi', 'thisFileNameHUFiles', 'indOfDash', 'convention')

    intCovMtxEpiAll = intCovMtxEpiAll + intCovMtxEpi;
    numerEpiAll = numerEpiAll + numerEpi;
end

denomEpiAll = (intCovMtxEpiAll + regMtxEpi);
selEstEpi_rep = denomEpiAll\numerEpiAll;

% find accessibility



if(size(intCovMtxEpiAll, 1) < 50000)
   [selcSitesEpi_rep] = findIndColsOfHmat(intCovMtxEpiAll);

   if(exist([dirNameAnalysisRep 'Estimates' chosenSlash], 'dir') == 0)
        mkdir([dirNameAnalysisRep 'Estimates' chosenSlash])
   end
   dlmwrite([dirNameAnalysisRep 'Estimates' chosenSlash  'AccessibilityMPLEpi_rep.txt'], selcSitesEpi_rep);
else
   disp('IntCovMtEpi size > 5000...case not handled, Accessibility not computed.')
end


% save results
fileNameSelEstEpi_rep = ['SelEstEpi_rep_gamma' num2str(priorConstSC) '_' num2str(priorConstEpi) '.txt'];

dlmwrite([dirNameAnalysisRep 'Estimates' chosenSlash fileNameSelEstEpi_rep], full(selEstEpi_rep))

