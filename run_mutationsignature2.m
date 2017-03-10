%identifying the number of mutational processes operative in a set of mutational catalogues
addpath('source/');
addpath('plotting/');
clc;

%% Open matlabpool
if ( matlabpool('size') == 0 )
    %matlabpool open; % opens the default matlabpool, if it is not already opened
    matlabpool open
end

%% Define parameters
iterationsPerCore = 100;
minNumberOfSignature = 1;
%maxNumberOfSignature = 15;
maxNumberOfSignature = 4;
stability = zeros(maxNumberOfSignature, 1);
reconstructionError = zeros(maxNumberOfSignature, 1);
disp(inputFile);
disp(allOutputFile);
%%inputFile = 'input/combineESCCHENAN_110samples_substitutions.mat';
%%allOutputFile = 'output/res_combineESCCHENAN_example2_110samples_substitutions.mat';

%% Sequentially deciphering signatures between minNumberOfSignature and maxNumberOfSignature
for totalSignatures = minNumberOfSignature : maxNumberOfSignature

    % Decipher the signatures of mutational processes from catalogues of mutations
    [input allProcesses allExposures idx processes exposures processStab processStabAvg] = ...
        decipherMutationalProcesses(iterationsPerCore, totalSignatures, inputFile, ...
            ['output/' pref '_' num2str(totalSignatures) '_signatures.mat'] );

    % Record the stability and average Frobenius reconstruction error
    stability(totalSignatures-minNumberOfSignature+1) = mean(processStabAvg);
    reconstructionError(totalSignatures-minNumberOfSignature+1) = norm(input.originalGenomes - processes*exposures, 'fro');

end

%% Plotting the stability and average Frobenius reconstruction error
try %% Some versions of MATLAB plotyy has a bug under linux with -nodisplay -nosplash -nodesktop options
  plotSignatureStabilityAndReconstruction(minNumberOfSignature:maxNumberOfSignature, stability, reconstructionError, input);
catch ME
  %% Do not do anything - just ignore the plot in order to save the final output daya
end
%% Saving the data
save(allOutputFile);
