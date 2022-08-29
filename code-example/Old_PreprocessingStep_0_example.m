 
% prepares  dirNames_X_Y.txt files needed in PreProcessingStep1

% usage: data folder should have the following structure
%        Data dir:  /.../dir1/[Patient_ID]/[Protein_name]/[FileNames].fasta
%        Analysis Dir should be different from data directiry


clc
clear all
close all

%-------------------- USER CONTROLLED INITIALIZATION ----------------------
FLAG_SaveFile = 1;


fileNameContainingDirPath = 'Set_dirNames_MPL_SimDataEpi.txt';
[dirNameData, dirNameAnalysis, dirNameResults] = setDirNamesMPLPipeline(fileNameContainingDirPath);

% These variables need to be set for given dataset
dataDirNameMain = dirNameData; %[dirNameData 'Set' num2str(thisSet) '\ng' num2str(ng) '_dT' num2str(dTStep) '_Tused' num2str(Tused) '_initStr' num2str(numStrainsInInitialPop) '\'];
analysisDirNameMain = dirNameAnalysis; %[dirNameAnalysis 'Set' num2str(thisSet) '\ng' num2str(ng) '_dT' num2str(dTStep) '_Tused' num2str(Tused) '_initStr' num2str(numStrainsInInitialPop) '\'];
resultsDirNameMain = dirNameResults; % [dirNameResults 'Set' num2str(thisSet) '\ng' num2str(ng) '_dT' num2str(dTStep) '_Tused' num2str(Tused) '_initStr' num2str(numStrainsInInitialPop) '\'];

str1 = 'dirNamesExample';% ['dirNamesSet' num2str(thisSet) '_ng' num2str(ng) '_dT' num2str(dTStep) '_Tused' num2str(Tused) '_initStr' num2str(numStrainsInInitialPop) '_'];
%--------------------------------------------------------------------------

% ------------------------- AUTO INITIALIZATION ---------------------------
% NO USER INPUT REQUIRED
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
if(exist(dirNameStr1Files, 'dir') == 0)
    mkdir(dirNameStr1Files)
end

dirNameDataThisPat = dataDirNameMain;
dirNameAnalysisThisPat = analysisDirNameMain; 
dirNameResultsThisPat = resultsDirNameMain; 

if(exist(dirNameAnalysisThisPat, 'dir') == 0)
    mkdir(dirNameAnalysisThisPat)
end

fileNameContainingDirPath = [str1 '.txt'];

if(FLAG_SaveFile == 1)
    fprintf('Saving .txt file containing names of data and analysis folders...')

    if(exist([dirNameStr1Files fileNameContainingDirPath], 'file') == 2)
        delete([dirNameStr1Files fileNameContainingDirPath])
    end
    fileID = fopen([dirNameStr1Files fileNameContainingDirPath],'w');
    for f = 1:3
       if(f == 1)
          fprintf(fileID,'%s\n', '% This file specifies directory paths to Data, Analysis');
          fprintf(fileID,'%s\n', '% the syntax is /dir1/dir2/dir3/ ');
          fprintf(fileID,'%s\n', '% for windows based system, the code will automatically reformat the path');
          fprintf(fileID,'%s\n', '%');
          fprintf(fileID,'%s\n', '% -------------------------------------------------------------------------');
       end
       if(f == 1)
           fprintf(fileID,'%s\n',['dirNameData=' dirNameDataThisPat]);
       elseif(f == 2)
           fprintf(fileID,'%s\n',['dirNameAnalysis=' dirNameAnalysisThisPat]);
       elseif(f == 3)
           fprintf(fileID,'%s\n',['dirNameResults=' dirNameResultsThisPat]);
       end
    end
    fclose(fileID);
    disp('done.')
elseif(FLAG_SaveFile == 0)
    disp('Warning: .txt data file not saved as FLAG_SaveFile flag not set.')
else
    disp('Error: case undefined')
    pause
end

