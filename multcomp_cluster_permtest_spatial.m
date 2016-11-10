function [Results] = multcomp_cluster_permtest_spatial(cond1_data, cond2_data, neighbourhood_matrix_filepath, varargin)
%
% This script receives paired-samples data for each condition and outputs corrected p-values and
% hypothesis test results based on a maximum cluster statistic permutation test.
% The permutation test in this script is based on the t-statistic, 
% but could be adapted to use with other statistics such as the trimmed mean.
% Clustering is performed across channels (spatial clustering).
%
% For details on this multiple testing procedure see:
% Bullmore, E. T., Suckling, J., Overmeyer, S., Rabe-Hesketh, S., 
% Taylor, E., & Brammer, M. J. (1999). Global, voxel, and cluster tests, 
% by theory and permutation, for a difference between two groups of 
% structural MR images of the brain. IEEE Transactions on Medical Imaging,
% 18, 32-42. doi 10.1109/42.750253
%
% Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing
% of EEG- and MEG-data. Journal of Neuroscience Methods, 164, 177-190.
% doi 10.1016/j.jneumeth.2007.03.024
%
% This function implements a conservative correction for p-values when
% using permutation tests, as described by Phipson & Smyth (2010).
% Permutation p-values should never be zero: Calculating exact p-values
% when permutations are randomly drawn. Statistical Applications in
% Genetics and Molecular Biology, 9, 39. doi 10.2202/1544-6115.1585
%
% Inputs:
%
%   cond1_data                      data from condition 1, 
%                                   a subjects x channels matrix
%
%   cond2_data                      data from condition 2, 
%                                   a subjects x channels matrix
%
%   neighbourhood_matrix_filepath   file path to the matlab file containing the
%                                   channel neighbourhood matrix channeighbstructmat
%
% Optional Inputs:
%
%   alpha                           uncorrected alpha level, default 0.05
%
%   iterations                      number of permutation samples to draw. 
%                                   At least 1000 is recommended for the 
%                                   p = 0.05 alpha level, and at least 5000 is
%                                   recommended for the p = 0.01 alpha level. 
%                                   This is due to extreme events at the tails
%                                   being very rare, needing many random permutations
%                                   to find enough of them.
%
%   clusteringalpha                 the significance threshold used to define
%                                   individual points within a cluster. 
%                                   Setting this to larger values (e.g.
%                                   0.05) will detect broadly distributed clusters,
%                                   whereas setting it to 0.01 will help detect 
%                                   smaller clusters that exhibit strong effects.
%
%   min_cluster_size                The minimum size of clusters which are
%                                   declared statistically significant.
%                                   This constraint does not apply to the
%                                   randomly-partitioned datasets that make
%                                   up the null permutation distribution.
%
%   yuen_t                          Set whether to use Yuen's t test
%                                   (robust version of Student's t). This  
%                                   1 = use Yuen's t, 0 = use Student's t
%
%   trimming                        Trimming parameter for Yuen's t.
%                                   Default is 20 (20% trimming)
%
% Outputs:
%
%   "Results" structure containing:
%
%   corrected_h                     vector of hypothesis tests for which statistical 
%                                   significance is defined by values above a threshold
%                                   of the (alpha_level * 100)th percentile
%                                   of the maximum statistic distribution.
%                                   1 = statistically significant, 0 = not statistically significant)
% 
%   t_values                        vector of t-values from the paired-samples 
%                                   tests at each channel.
%
%   cluster_channel_indices         cell array of channel indices in each cluster.
%                                   For easy separating of clusters when plotting 
%                                   and interpreting results.
% 
%   cluster_mass_vector             vector of summed t-values within each cluster,
%                                   i.e. the mass of each cluster. Used for comparing
%                                   against the null cutoff for NHST.
%
%   cluster_mass_null_cufoff_pos    the (1-alpha/2)th percentile of maximum cluster
%                                   masses returned by the cluster-based permutation
%                                   test. i.e. the null threshold for positive direction 
%                                   effects.
%
%   cluster_mass_null_cufoff_pos    the (alpha/2)th percentile of maximum cluster
%                                   masses returned by the cluster-based permutation
%                                   test. i.e. the null threshold for negative direction 
%                                   effects.
%
%   cluster_p                       p-value of cluster, calculated as the
%                                   number of permutation sample maximum
%                                   cluster masses larger than the observed
%                                   cluster mass + 1 divided by the total
%                                   number of permutations + 1. This is the
%                                   conservative p-value correction in
%                                   Phipson & Smyth (2010).
%
% Example:      [RESULTS] = multcomp_cluster_permtest_spatial(cond1_data, cond2_data, 'expected_chanlocs.mat', 'alpha', 0.05, 'iterations', 1000, 'clusteringalpha', 0.05, 'yuen_t', 1)
%
%
% Copyright (c) 2016 Daniel Feuerriegel
%
%
% This file is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
%
% NOTE: In this script a 'step' refers to a channel (i.e. a separate
% analysis). This is carried over from the temporal clustering script which
% used 'step' to mean a time-window of interest.
%

%% Handling variadic inputs
% Define defaults at the beginning
options = struct(...
    'alpha', 0.05,...
    'iterations', 1000,...
    'clusteringalpha', 0.05,...
    'min_cluster_size', 1,...
    'yuen_t', 0,...
    'trimming', 20);
    

% Read the acceptable names
option_names = fieldnames(options);

% Count arguments
n_args = length(varargin);
if round(n_args/2) ~= n_args/2
   error([mfilename ' needs property name/property value pairs'])
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inp_name = lower(pair{1}); % make case insensitive

   % Overwrite default options
   if any(strcmp(inp_name, option_names))
      options.(inp_name) = pair{2};
   else
      error('%s is not a recognized parameter name', inp_name)
   end
end
clear pair
clear inp_name

% Renaming variables for use below:
alpha_level = options.alpha;
n_iterations = options.iterations;
clustering_alpha = options.clusteringalpha;
min_cluster_size = options.min_cluster_size;
use_yuen = options.yuen_t;
trimming = options.trimming;
clear options;

%% Seed the random number generator
rng('shuffle'); % Seed random number generator based on computer clock

%% Loading Channel Locations
load(neighbourhood_matrix_filepath); % Loads channeighbstructmat matrix of channel neighbourhoods
% This matrix marks which electrodes are adjacent/connected to which others. 
% (1 = connected/ 0 = not connected)


%% Cluster-based multiple testing correction

% Checking whether the number of steps of the first and second datasets are equal
if size(cond1_data, 2) ~= size(cond2_data, 2)
   error('Condition 1 and 2 datasets do not contain the same number of channels!');
end
if size(cond1_data, 1) ~= size(cond2_data, 1)
   error('Condition 1 and 2 datasets do not contain the same number of subjects!');
end

% Generate difference scores between conditions
diff_scores = cond1_data - cond2_data;

n_subjects = size(diff_scores, 1); % Calculate number of subjects
n_total_comparisons = size(diff_scores, 2); % Calculating the number of comparisons

% Perform t-tests at each step
uncorrected_h = zeros(1, n_total_comparisons); % preallocate
uncorrected_t = zeros(1, n_total_comparisons); % preallocate
uncorrected_t_sign = zeros(1, n_total_comparisons); % preallocate

if use_yuen == 0
    [uncorrected_h, ~, ~, extra_stats] = ttest(diff_scores, 0, 'Alpha', clustering_alpha);
    uncorrected_t = extra_stats.tstat; % Recording t statistic for each test
elseif use_yuen == 1
    for step = 1:n_total_comparisons
        [uncorrected_h(step), ~, ~, t_yuen, ~, ~, ~, ~] = yuend_ttest(cond1_data(:, step), cond2_data(:, step), trimming, clustering_alpha);
        uncorrected_t(step) = t_yuen; % Recording t statistic for each test
    end
end

for step = 1:n_total_comparisons
    % Marking the direction of effects
    if uncorrected_t(step) < 0
        uncorrected_t_sign(step) = -1;
    else
        uncorrected_t_sign(step) = 1;  
    end % of if uncorrected_t(step)
end % of for step

% Generate maximum cluster mass distribution from randomly-partitioned datasets
t_stat = zeros(n_total_comparisons, n_iterations); % Preallocate
max_pos_cluster_mass = zeros(1, n_iterations); % Preallocate
max_neg_cluster_mass = zeros(1, n_iterations); % Preallocate
cluster_perm_test_h = zeros(n_total_comparisons, n_iterations); % Preallocate
t_sign = zeros(n_total_comparisons, n_iterations); % Preallocate

for iteration = 1:n_iterations
    % Draw a random partition sample for each test
    temp_dataset_1 = zeros(n_subjects, n_total_comparisons);
    temp_dataset_2 = zeros(n_subjects, n_total_comparisons);
    temp_diffscores = zeros(n_subjects, n_total_comparisons);
    
    for step = 1:n_total_comparisons 
        % Randomly switch the sign of difference scores (equivalent to
        % switching labels of conditions/random partitioning)
        temp_signs = (rand(n_subjects, 1) > .5) * 2 - 1; % Switches signs of labels
        
        temp_dataset_1(temp_signs == 1, step) = cond1_data(temp_signs == 1, step);
        temp_dataset_1(temp_signs == -1, step) = cond2_data(temp_signs == -1, step);
        temp_dataset_2(temp_signs == 1, step) = cond2_data(temp_signs == 1, step);
        temp_dataset_2(temp_signs == -1, step) = cond1_data(temp_signs == -1, step);
        
        temp_diffscores(1:n_subjects, step) = temp_dataset_1(1:n_subjects, step) - temp_dataset_2(1:n_subjects, step);
    end % of for step
        
    % Perform t tests (Student's or Yuen's paired-samples t)
    if use_yuen == 0
        [cluster_perm_test_h(:, iteration), ~, ~, temp_stats] = ttest(temp_diffscores, 0, 'Alpha', clustering_alpha);
        t_stat(:, iteration) = temp_stats.tstat; % Get t statistic
    elseif use_yuen == 1
        for step = 1:n_total_comparisons
            [cluster_perm_test_h(step, iteration), ~, ~, t_yuen, ~, ~, ~, ~] = yuend_ttest(temp_dataset_1(:, step), temp_dataset_2(:, step), trimming, clustering_alpha);
            t_stat(step, iteration) = t_yuen; % Recording t statistic for each test
        end
    end
    
    
        % Marking the sign of each t statistic to avoid clustering pos
        % and neg significant results
        for step = 1:n_total_comparisons
            if t_stat(step, iteration) < 0;
               t_sign(step, iteration) = -1; 
            else
               t_sign(step, iteration) = 1; 
            end % of if t_stat statement
        end
    % Initialise a matrix recording whether a channel has been assessed
    channel_processed = zeros(1, n_total_comparisons);
    
    % Initialise a cell array to record which channels belong to which cluster
    cluster_channel_indices = {};
    
    % Identify clusters and generate a maximum cluster statistic
    cluster_mass_vector = [0]; % Resets vector of cluster masses
    cluster_counter = 0;
    
    for step = 1:n_total_comparisons    
        
        % Check whether channel is statistically significant and if it has
        % already been processed
        if cluster_perm_test_h(step, iteration) == 1 && channel_processed(step) == 0
        
            % Mark that this channel has been processed
            channel_processed(step) = 1;
        
            % Increase the cluster counter index
            cluster_counter = cluster_counter + 1;
            
            % Add channel to vector of cluster channels
            cluster_channel_indices{cluster_counter}(1) = step;
            
            % Mark neighboring channels
            neighbouring_channels = []; % Reset neighbouring channels vector
            channel_counter = 1; % Reset counter for neighbouring channels
            for i = 1:n_total_comparisons
               % If channel is a neighbour, not yet processed and
               % statistically significant in the same direction as the
               % other channels in the cluster
               if channeighbstructmat(step, i) == 1 && cluster_perm_test_h(i, iteration) == 1 && channel_processed(i) == 0 && t_sign(step, iteration) == t_sign(i, iteration)

                   % Mark that this channel has been processed
                   channel_processed(i) = 1;

                   % Add channel to vector of cluster channels
                   cluster_channel_indices{cluster_counter}(end + 1) = i;

                   % Mark it as a neighbouring channel for searching below
                   neighbouring_channels(channel_counter) = i;
                   channel_counter = channel_counter + 1;
                   
               end % of if channneighbstructmat
            end % of for i = 1:nChannels
                   
            % Scan for neighboring channels to add to that cluster
            keep_searching = 1; % Initialise variable for while loop

            while keep_searching == 1

                keep_searching = 0;

                % Set channels to search
                chans_to_search = neighbouring_channels;
                neighbouring_channels = []; % clear out neighbouring channels vector
                neighbouring_channel_counter = 1;

                for chanind = chans_to_search
                    for k = 1:n_total_comparisons
                        % Check whether there are neighbouring channels with sig.
                        % results in the same direction
                        if cluster_perm_test_h(k, 1) == 1 && channel_processed(k) == 0 && channeighbstructmat(chanind, k) == 1 && t_sign(step, iteration) == t_sign(k, iteration)

                            % Set to keep searching for other channels in next loop iteration
                            keep_searching = 1; 

                            % Mark channel as having been already processed
                            channel_processed(k) = 1;

                            % Add any sig neighbouring channels to the channel vector
                            % for the cluster 
                            cluster_channel_indices{cluster_counter}(end + 1) = k;

                            % Set to search these channels for further
                            % neighbours in the next iteration of the while loop
                            neighbouring_channels(neighbouring_channel_counter) = k;
                            neighbouring_channel_counter = neighbouring_channel_counter + 1;

                        end % of if cluster_perm_test_h
                    end % of for k
                end % of for chanind
            end % of while keepsearching

            % Sum the t-values of the vector of neighbouring sig. channels to
            % make the cluster statistic. Record in a matrix
            cluster_mass_vector(cluster_counter) = sum(t_stat(cluster_channel_indices{cluster_counter})); 
            
        else % If the channel was not statistically significant with respect to clustering alpha
        
            % Mark channel as having been already processed
            channel_processed(step) = 1;
            
        end % of if cluster_perm_test_h statement
        
    end % of for step
    
    % Find the maximum cluster mass for positie and negative values
    max_pos_cluster_mass(iteration) = max(cluster_mass_vector);
    max_neg_cluster_mass(iteration) = min(cluster_mass_vector);
   
end % of iteration loop

% Calculate the 95th percentiles of maximum cluster mass values (used as decision
% critieria for statistical significance)
cluster_mass_null_cutoff_pos = prctile(max_pos_cluster_mass, ((1 - alpha_level /2) * 100));
cluster_mass_null_cutoff_neg = prctile(max_neg_cluster_mass, ((alpha_level /2) * 100));

%% Calculate cluster masses in the actual (non-permutation) tests

% Initialise a matrix recording whether a channel has been assessed
channel_processed = zeros(1, n_total_comparisons);

% Initialise a cell array to record which channels belong to which cluster
cluster_channel_indices = {};

% Initialise cluster variables
cluster_mass_vector = [0]; % Resets vector of cluster masses
cluster_counter = 0;

% Initialise vector of cluster-corrected significant channels
cluster_corrected_sig_steps = zeros(1, n_total_comparisons);

for step = 1:n_total_comparisons    
        
    % Check whether channel is statistically significant and if it has
    % already been processed
    if uncorrected_h(step) == 1 && channel_processed(step) == 0

        % Mark that this channel has been processed
        channel_processed(step) = 1;

        % Increase the cluster counter index
        cluster_counter = cluster_counter + 1;

        % Add channel to vector of cluster channels
        cluster_channel_indices{cluster_counter}(1) = step;

        % Mark neighboring channels
        neighbouring_channels = []; % Reset neighbouring channels vector
        channel_counter = 1; % Reset counter for neighbouring channels
        for i = 1:n_total_comparisons % cycle through all channels
           % If channel is a neighbour, not yet processed and
           % statistically significant in the same direction as the
           % other channels in the cluster
           if channeighbstructmat(step, i) == 1 && uncorrected_h(i) == 1 && channel_processed(i) == 0 && uncorrected_t_sign(step) == uncorrected_t_sign(i)

               % Mark that this channel has been processed
               channel_processed(i) = 1;

               % Add channel to vector of cluster channels
               cluster_channel_indices{cluster_counter}(end + 1) = i;

               % Mark it as a neighbouring channel for searching below
               neighbouring_channels(channel_counter) = i;
               channel_counter = channel_counter + 1;

           end % of if channneighbstructmat
        end % of for i = 1:nChannels

        % Scan for neighboring channels to add to that cluster
        keep_searching = 1; % Initialise variable for while loop

        while keep_searching == 1
            
            keep_searching = 0;

            % Set channels to search
            chans_to_search = neighbouring_channels;
            neighbouring_channels = []; % clear out neighbouring channels vector
            neighbouring_channel_counter = 1;

            for chanind = chans_to_search
                for k = 1:n_total_comparisons
                    % Check whether there are neighbouring channels with sig.
                    % results in the same direction
                    if uncorrected_h(k) == 1 && channel_processed(k) == 0 && channeighbstructmat(chanind, k) == 1 && uncorrected_t_sign(step) == uncorrected_t_sign(k)

                        % Set to keep searching for other channels in next loop iteration
                        keep_searching = 1; 

                        % Mark channel as having been already processed
                        channel_processed(k) = 1;

                        % Add any sig neighbouring channels to the channel vector
                        % for the cluster 
                        cluster_channel_indices{cluster_counter}(end + 1) = k;

                        % Set to search these channels for further
                        % neighbours in the next iteration of the while loop
                        neighbouring_channels(neighbouring_channel_counter) = k;
                        neighbouring_channel_counter = neighbouring_channel_counter + 1;

                    end % of if cluster_perm_test_h
                end % of for k
            end % of for chanind
        end % of while keepsearching

        % Sum the t-values of the vector of neighbouring sig. channels to
        % make the cluster statistic. Record in a matrix
        cluster_mass_vector(cluster_counter) = sum(uncorrected_t(cluster_channel_indices{cluster_counter})); 

    else % If the channel was not statistically significant
        
        % Mark channel as having been already processed
        channel_processed(step) = 1;
        
    end % of if cluster_perm_test_h statement
end % of for step

% Determine whether the cluster masses are above the null cutoff and meet
% the minimum cluster size constraint
cluster_sig = zeros(length(cluster_mass_vector), 1); % Preallocate

for cluster_no = 1:length(cluster_mass_vector);
    if cluster_mass_vector(cluster_no) > cluster_mass_null_cutoff_pos || cluster_mass_vector(cluster_no) < cluster_mass_null_cutoff_neg
        if length(cluster_channel_indices{cluster_no}) >= min_cluster_size
            cluster_sig(cluster_no) = 1; % mark cluster as statistically significant
            cluster_corrected_sig_steps(cluster_channel_indices{1, cluster_no(:)}) = 1;
        end 
    end % of if cluster_mass_vector
end % of for cluster_no

% Create vectors of multiple comparisons corrected hypothesis tests and the
% t-values derived from tests at each channel
corrected_h = cluster_corrected_sig_steps;
t_values = uncorrected_t;

% Calculating a p-value for each cluster
% p-values are corrected according to Phipson and Smyth (2010) methods
cluster_p = ones(length(cluster_mass_vector), 1); % Preallocate p-values

for cluster_no = 1:length(cluster_mass_vector)
    
    % Calculate the number of permutation samples with cluster masses
    % larger than the observed cluster mass for a given cluster
    if cluster_mass_vector(cluster_no) > 0 % If a positive difference
        b = sum(abs(max_pos_cluster_mass) >= abs(cluster_mass_vector(cluster_no)));
    elseif cluster_mass_vector(cluster_no) < 0 % If a negative difference
        b = sum(abs(max_neg_cluster_mass) >= abs(cluster_mass_vector(cluster_no)));
    end    
    p_t = (b + 1) / (n_iterations + 1); % Calculate conservative version of p-value as in Phipson & Smyth, 2010
    
    cluster_p(cluster_no) = p_t * 2; % Multiplied by 2 as two-tailed
       
end % of for cluster_no


% 
% for cluster_no = 1:length(cluster_mass_vector);
%     if cluster_mass_vector(cluster_no) > 0
%         cluster_p(cluster_no) = mean(max_pos_cluster_mass(:) >= cluster_mass_vector(cluster_no)) * 2; % Multiply by 2 for two-tailed
%     elseif cluster_mass_vector(cluster_no) < 0
%         cluster_p(cluster_no) = mean(max_neg_cluster_mass(:) <= cluster_mass_vector(cluster_no)) * 2; % Multiply by 2 for two-tailed
%     end
% end % of for cluster_no

% Copying the results into a structure to output
Results.uncorrected_h = uncorrected_h;
Results.corrected_h = corrected_h;
Results.t_values = t_values;
Results.cluster_sig = cluster_sig;
Results.cluster_channel_indices = cluster_channel_indices;
Results.cluster_mass_vector = cluster_mass_vector;
Results.cluster_mass_null_cutoff_pos = cluster_mass_null_cutoff_pos;
Results.cluster_mass_null_cutoff_neg = cluster_mass_null_cutoff_neg;
Results.cluster_p = cluster_p;
end % of function
