 
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


% specify complete paths to the folder that will contain Data, Analysis, and Results
dataDirNameMain = 'D:/epistasis-inference-example/Example_Data_1/';
analysisDirNameMain = 'D:/epistasis-inference-example/Example_Analysis_1/';
resultsDirNameMain = 'D:/epistasis-inference-example/Example_Results_1/';

% this is the filename where the above paths will be stored. This file is
% located at ...\epistasis-inference-example\Data_Misc\dirNameFiles\
% This filename should be unique for each dataset being analyzed
fileNameContainingDirPathCell{1} = 'dirNamesExample';


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
analysisMPL_Epi_v2_outerFunction(priorConstSC, priorConstEpi, recombProb, refernceSequence, ...
    dataDirNameMain, analysisDirNameMain, resultsDirNameMain, fileNameContainingDirPathCell, ...
    fileNameContainingMutProb, Lin)


% estimated selection coefficients are in
% ...\epistasis-inference-example\Example_Analysis\Estimates\SelEstEpi_MC1_synth_Ito_gammaX_Y.txt
% where X and Y are the user provided numerical values of priorConstSC and
% priorConstEpi 

% FORMAT of SelEstEpi_MC1_synth_Ito_gammaX_Y.txt: 
% It has Lin*(Lin+1)/2 entries where Lin is the user provided length of the
% aligned sequence

% The first Lin entries are the selection coefficients (s_i) while the rest
% are pairwise epistasis terms (s_ij) in the following order:

% s_1 s_2 ... s_Lin s_12 s_13 ... s_1Lin s_23 s_24 ... s_2Lin ... s_(Lin-1)Lin

%% 
plot_MPL_outerFunction(Lin, {dataDirNameMain}, {analysisDirNameMain}, priorConstSC, priorConstEpi)