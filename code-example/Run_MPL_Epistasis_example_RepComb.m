 
% prepares  dirNames_X_Y.txt files needed in PreProcessingStep1

% usage: data folder should have the following structure
%        Data dir:  /.../dir1/[Patient_ID]/[Protein_name]/[FileNames].fasta
%        Analysis Dir should be different from data directiry


clc
clear all
close all

%-------------------- USER CONTROLLED INITIALIZATION ----------------------

priorConstSC = 1; % this is the strength of the selection coefficient regularization parameter
priorConstEpi = 1; % this is the strength of the epistasis terms regularization parameter

numRep = 5; % number of replicates to combine

% specify complete paths to the folder that will contain Data, Analysis, and Results
for i = 1:numRep
    dataDirNameMainCell{i} = ['D:/epistasis-inference-example/Example_Data_' num2str(i) '/'];
    analysisDirNameMainCell{i} = ['D:/epistasis-inference-example/Example_Analysis_' num2str(i) '/'];
    resultsDirNameMainCell{i} = ['D:/epistasis-inference-example/Example_Results_' num2str(i) '/'];


    % this is the filename where the above paths will be stored. This file is
    % located at ...\epistasis-inference-example\Data_Misc\dirNameFiles\
    % This filename should be unique for each dataset being analyzed
    fileNameContainingDirPathCell{i} = ['dirNamesRepCombExample_1_' num2str(i)];
end
dirNameAnalysisRep = ['D:/epistasis-inference-example/Example_Analysis_RepComb_1_5/'];

Lin = 5; % specify the length of the aligned sequence

% specify NT reference sequence
refernceSequence = repmat('A', 1, Lin);


% this is the filename which contains the NT-to-NT mutation probability. It
% must be located in the folder
% ...\epistasis-inference-example\Data_Misc\MutationProbabilities\ 
fileNameContainingMutProb = 'MutProb_SyntheticData_1eminus4.txt';

recombProb = 1e-4; % recombination probability


%--------------------------------------------------------------------------
%------------ NO USER INPUT REQUIRED ------------


%%

% Analysis
analysisMPL_Epi_v2_rep_outerFunction(priorConstSC, priorConstEpi, recombProb, refernceSequence, ...
    dataDirNameMainCell, analysisDirNameMainCell, resultsDirNameMainCell, fileNameContainingDirPathCell, ...
    fileNameContainingMutProb, Lin, numRep, dirNameAnalysisRep)




%%
plot_MPL_rep_outerFunction(Lin, numRep, dataDirNameMainCell, analysisDirNameMainCell, dirNameAnalysisRep, priorConstSC, priorConstEpi)
