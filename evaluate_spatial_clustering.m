% evaluate_spatial_clustering.m
%
% Used to evaluate the weak FWER of the cluster-based script. A sanity
% check to make sure the function provides adequate FWER control
% 
% Written by DF on 10/16

% Housekeeping
% clear all;
close all;

% See random number generator based on computer clock
rng('shuffle');

% Simulation variables
useRealData = 0; % 1 = Use permuted data from Unexpected Oddball Exp 1; 0 = Use simulated data
nIterations = 100; % Number of iterations
nIterationsPermtest = 1000; % Number of iterations to use within multiple comparisons correction function (default = 1000)
useRobust = 0; % 1 = Use yuen's t test; 0 = Use standard t test
alphaLevel = 0.05; % Testwise alpha level
clusteringAlpha = 0.05; % Alpha level used to determeine whether individual electrodes will be included within a cluster

% Determine minimum cluster size (minimum number of channels in cluster to be declared a sig. effect)
minClusterSizes = [1:10]; 

% Choose sample size and number of channels
sampleSize = 22;
nChannels = 128; % Number of channels - 128 used for testing biosemi data with default chanlocs file

type1Errors = zeros(nIterations, length(minClusterSizes)); % Reset number of type 1 errors

% Load datasets to use for simulations
if useRealData == 1
    load('sumHarm_19');
    realEEGDataset1 = sumHarm_19(1:nChannels, 1:sampleSize);
    realEEGDataset1 = permute(realEEGDataset1, [2 1]);

    load('sumHarm_29');
    realEEGDataset2 = sumHarm_29(1:nChannels, 1:sampleSize);
    realEEGDataset2 = permute(realEEGDataset2, [2 1]);
end

% Evaluate the FWER
for iteration = 1:nIterations

    fprintf(['running iteration ' int2str(iteration) '...\n']);
    
    if useRealData == 0 % Randomly generate datasets
        % Randomly generate data from conditions 1 and 2
        dataset_1 = zeros(sampleSize, nChannels);
        dataset_2 = zeros(sampleSize, nChannels);
        for i = 1:nChannels
            dataset_1(:,i) = randn(sampleSize, 1);
            dataset_2(:,i) = randn(sampleSize, 1);
        end
    elseif useRealData == 1 % Randomly permute 2 real EEG datasets
        % Preallocate matrices
        dataset_1 = zeros(sampleSize, nChannels);
        dataset_2 = zeros(sampleSize, nChannels);
        
        for channel = 1:nChannels  
            temp_signs = (rand(sampleSize, 1) > .5) * 2 - 1; % Switches signs of labels
            dataset_1(temp_signs == 1, channel) = realEEGDataset1(temp_signs == 1, channel);
            dataset_1(temp_signs == -1, channel) = realEEGDataset2(temp_signs == -1, channel);
            dataset_2(temp_signs == 1, channel) = realEEGDataset2(temp_signs == 1, channel);
            dataset_2(temp_signs == -1, channel) = realEEGDataset1(temp_signs == -1, channel);
        end % of for channel 
    end % of if useRealData
        
    % Perform the mass-univariate testing (robust version)
    if useRobust == 1
        
        [Results] = multcomp_cluster_permtest_spatial(dataset_1, dataset_2, 'expected_chanlocs.mat', 'alpha', alphaLevel, 'clusteringalpha', clusteringAlpha, 'iterations', nIterationsPermtest, 'yuen_t', 1);

    elseif useRobust == 0
        
        [Results] = multcomp_cluster_permtest_spatial(dataset_1, dataset_2, 'expected_chanlocs.mat', 'alpha', alphaLevel, 'clusteringalpha', clusteringAlpha, 'iterations', nIterationsPermtest);
    
    end
    
    % Check for sig. clusters and cluster size constraints. Declare if a
    % Type 1 error has occured (weak FWER error control)
    for clusterSizeThreshold = 1:length(minClusterSizes)
        for j = 1:length(Results.cluster_mass_vector)
            if Results.cluster_mass_vector(j) > Results.cluster_mass_null_cutoff_pos && length(Results.cluster_channel_indices{j}) >= minClusterSizes(clusterSizeThreshold)
               type1Errors(iteration, clusterSizeThreshold) = 1;              
            end  % of if statement
            
            if Results.cluster_mass_vector(j) < Results.cluster_mass_null_cutoff_neg && length(Results.cluster_channel_indices{j}) >= minClusterSizes(clusterSizeThreshold)
               type1Errors(iteration, clusterSizeThreshold) = 1;              
            end  % of if statement
        end % for j = 1:length(cluster_mass_vector)
    end % of for clusterSizeThreshold
end % of for iteration

for i = 1:length(minClusterSizes)
    totalType1Errors(i) = sum(type1Errors(:,i)); % Calculate total number of false positives)
    weakFWER(i) = totalType1Errors(i) / nIterations; % Number of type 1 errors divided by number of iterations
end

% Create save directory if doesn't exist
if ~exist('Workspace Saves', 'dir') 
    mkdir('Workspace Saves');
end

% Save the workspace to a .mat file

save(['Workspace Saves/realdata_' int2str(useRealData) '_samplesize_' int2str(sampleSize) '_iterations_' int2str(nIterations) '_' datestr(now, 30)]); 


