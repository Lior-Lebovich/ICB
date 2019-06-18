%% Inhibitory selectivity - population analysis:

% The columns of the file below correspond to:
% 1. nu_E_U, % 2. nu_E_D, % 3. nu_I_U, % 4. nu_I_D,
% where nu_Type_alpha is the population-average firing rates at the time of
% the decision. 
% The rows are sorted by network (1..10) --> trials (1..500).

addpath('analyses'); addpath('data');

fileName = 'selectivity-choice-test100-k=800-gee=0.3-gei=1.5-gie=2.0-gii=2-t=9s-ies=0.06-iis=0.04-ie=0.06-ii=0.04-ne=32000-ni=80.4-NTRIALS=500-NARCH=10-SELECT-NET01'
selectivity_mat = load(fileName);

% for each trial, 1 if nu_Type_U > nu_Type_D and -1 if nu_Type_U < nu_Type_D:
winning_E = 2 * ( selectivity_mat(:,1) > selectivity_mat(:,2) ) - 1;
winning_I = 2 * ( selectivity_mat(:,3) > selectivity_mat(:,4) ) - 1;

% compute the probability that the winning population has a higher 
% population-average firing rate, for both its excitatory and inhibitory
% neurons:
sum( winning_E == winning_I ) / length( winning_E )
