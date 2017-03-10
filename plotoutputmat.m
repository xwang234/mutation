%modified from example3.m

%% Example of reading a previously saved output file and plotting the
%% signatures of mutational processes and their respective exposures
%
% Ludmil B. Alexandrov
% Cancer Genome Project
% Wellcome Trust Sanger Institute
% la2@sanger.ac.uk
%
% This software and its documentation are copyright 2012 by the
% Wellcome Trust Sanger Institute/Genome Research Limited. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Wellcome Trust Sanger Institute nor Genome Research Limited 
% is responsible for its use, misuse, or functionality.

addpath('source/');
addpath('plotting/');
clc;

%% Define parameters
%outputFileName = 'output/res_example_1_21_WTSI_BRCA_whole_genome_substitutions.mat';
%outputFileName = 'output/res_henan_example1_6samples_substitutions.mat';
%outputFileName = 'output/res_henan_example2_6samples_substitutions_2_signatures.mat';
%outputFileName = 'output/res_otherdata_example1_104samples_substitutions.mat';
%outputFileName = 'output/escc_wgs_3_signatures.mat';
%outputFileName = 'output/escc_wes_3_signatures.mat';
%outputFileName = 'output/henan_varscan1_3_signatures.mat';
%outputFileName = 'output/dulak_varscan_3_signatures.mat';
%outputFileName = 'output/dulak_varscan2_3_signatures.mat';
%outputFileName = 'output/henan_varscan2_3_signatures.mat'
%outputFileName = 'output/escc_varscan2_3_signatures.mat'

if ( exist(outputFileName, 'file') == 0)
   disp(['File ' outputFileName ' cannot be found! Please run example1.m prior to running example3.m!']);
   return;
end

%% Loading the data
load(outputFileName);

%% Plotting the signatures of mutational processes
plotSignatures(processes, input, allProcesses, idx, processStabAvg);

%% Plotting the signature exposures
plotSignaturesExposureInSamples(exposures, input);
