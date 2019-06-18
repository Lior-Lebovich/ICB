%%  Identify the potential contribution of stochasticity at minutes’ time-
%   scale: Testing whether ICBs differ between the first and second halves 
%   of the experiments. 


%% Calculating the p.val matrix of all possible pairs of ICB in the 1st and
% 2nd half of the experiment ( ICB_half = 0:0.1:1).
% Running a permutation test on the order of trials to test:
% H0: ICB_1st-ICB_2nd = 0 vs. H1: ICB_1st-ICB_2nd <> 0 is estimated by

addpath('analyses','data');

nTrials = 20; % impossible trials
nPosValsHalf = (1 + 0.5 * nTrials);
binEdges = -1.05:0.1:0.05; % becuase possible |ICB_1st-ICB_2nd| = 0:0.1:1
nSim = 1e6;
pVals = zeros(nPosValsHalf,nPosValsHalf);
for i = 1:nPosValsHalf % p in 1st half = 0:0.1:1
    for j = 1:nPosValsHalf % p in 2nd half = 0:0.1:1
        % compute |ICB_1st - ICB_2nd| for each ICB_1st and ICB_2nd:
        realDiff = abs( (2/nTrials) * ( (i-1) - (j-1) ) );
        % permutation on |ICB_1st - ICB_2nd|:
        diff = zeros(1,nSim);
        for r = 1:nSim
            a = datasample( [ones(1,i+j-2) zeros(1,nTrials+2-i-j)], ...
                nTrials, 'Replace', false );
            diff(r) = abs( mean( a(1:(0.5*nTrials)) ) - ...
                mean( a((1+0.5*nTrials):nTrials) ) );
        end
        % compute the p.value:
        [l1,l2] = histcounts(-abs(diff), 'normalization',...
            'cdf', 'binEdges', binEdges);
        loc = find(l2<=-realDiff, 1, 'last' );
        if min(size(loc)) == 0
            pVals(i,j) = 0;
        elseif loc == length(l2)
            pVals(i,j) = 1;
        else
            pVals(i,j) = l1( loc );
        end
    end
end


%% Vertical bisection task:

nParticipantsBL = 100;

% load data:
allTrials = xlsread('BL_vertical_raw.xlsx', 1, 'C2:D12001');
impTrials = allTrials( allTrials(:,1) == 0 , 2 ); % consider only impossible trials

% compute, for each participant, pUp in the 1st and 2nd halves of the
% experiment:
p10_1st_2nd_20 = zeros(nParticipantsBL,3);
real_pVals = zeros(nParticipantsBL,1);
for i = 1:nParticipantsBL
    p10_1st_2nd_20(i,1) = mean(impTrials( ...
        (1 +(nTrials*(i-1))) : ((0.5*nTrials)+(nTrials*(i-1))) ) == 1); % 1st half
    p10_1st_2nd_20(i,2) = mean(impTrials( ...
        ((1+0.5*nTrials)+(nTrials*(i-1))) : (nTrials+(nTrials*(i-1))) ) == 1); % 2nd half
    p10_1st_2nd_20(i,3) = mean(impTrials( ...
        (1 +(nTrials*(i-1))) : (nTrials+(nTrials*(i-1))) ) == 1); % all trials
    % load relevant p.values:
    real_pVals(i) = pVals(  (0.5*nTrials) * p10_1st_2nd_20(i,1) + 1 ,...
        (0.5*nTrials) * p10_1st_2nd_20(i,2) + 1  );
end

% comute #participants for which pUp in 1st is significantly different than
% pUp in 2nd half of the experiment:
significant_BL = sum(real_pVals<=0.05);


%% Motor task:

nParticipantsMotor = 20;
nPairsMotor = 10;

% load data:
load('motor_raw');

% compute, for each participant-pair (20x10), pCW in the 1st and 2nd halves
% of the experiment:
p10_1st = zeros(1,nParticipantsMotor*nPairsMotor);
p10_2nd = zeros(1,nParticipantsMotor*nPairsMotor);
real_pVals = zeros(1,nParticipantsMotor*nPairsMotor);
for i = 1:nParticipantsMotor % running over participants
    p10_1st( (1+(0.5*nTrials)*(i-1)) : ((0.5*nTrials)*i) ) = ...
        mean( choices_all( 1:(0.5*nTrials),  1:(0.5*nTrials), i) ); % 1st half
    p10_2nd( (1+(0.5*nTrials)*(i-1)) : ((0.5*nTrials)*i) ) = ...
        mean( choices_all( (1+0.5*nTrials):nTrials, 1:(0.5*nTrials), i) ); % 2nd half
end

% load relevant p.values:
for i = 1:(nParticipantsMotor*nPairsMotor) % running over paricipantXpairs
    real_pVals(i) = pVals(  (0.5*nTrials)*p10_1st(i)+1 , ...
        (0.5*nTrials)*p10_2nd(i)+1  );
end

% comute #participants-pairs for which pCW in 1st is significantly
% different than pCW in 2nd half of the experiment:
significant_motor = sum(real_pVals<=0.05);


%% Compute percent significant:

percent_significant = (significant_BL + significant_motor) / ...
    (nParticipantsBL + nParticipantsMotor*nPairsMotor)