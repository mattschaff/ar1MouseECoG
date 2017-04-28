function [ avg_zscore, zscores ] = return_avg_zcore( true_parameters, trial_parameters )
%UNTITLED Summary of this function goes here
%   Inputs
    % true_parameters = (matrix) 1 X time 
    % trial_parameters = (matrix) trial_results  X time

% error if no parameters
if exist('true_parameters', 'var') == 0
    error('true_parameters is not set.');
end
if exist('trial_parameters', 'var') == 0
    error('trial_parameters is not set.');
end

%
%trial_parameters = [3 3 4 5 6 7; 2 3 4 5 6 7; 2 3 4 5 6 7; 2 3 5 5 6 7; 5 3 4 5 6 7;];
%true_parameters = [2 2 2 2 2 2];
trials_mean = nanmean(trial_parameters, 1);
trials_stdev = nanstd(trial_parameters);
zscores = abs((true_parameters - trials_mean)./(trials_stdev + .000001));
avg_zscore = nanmean(zscores); 
end

