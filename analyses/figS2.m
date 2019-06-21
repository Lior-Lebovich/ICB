%% supp. Fig. S2:

addpath('analyses','data');

% pPrev := the probability to repeat in a vertical impossible trial, the 
% choice (up or down) made in the previous vertical trials. Note that this
% corresponds to trial number k vs. k-4 (because an inpossible trial is
% always first in a triplet and it is preceded by 3 horizontal trials.

% The script below computes the pdf's, 90% confidence intervals and 100% 
% confidence intervals for each possible pair of pUp(impossible) and 
% pUp(before impossible);
% Relevant variables are saved to: prob4pPrevNotBiggerThanP2_data_SHORT.mat
prob4pPrevNotBiggerThanP2_short 


% Compute pUp in impossible trials and in preceding trials
nTrials = 120;
nImpTrials = 20;
data = xlsread('BL_vertical_raw.xlsx', 1, 'B2:E12001'); % columns: 1. trial 
% no. 2. deviation 3. answer (1=up, 0=down) 4. RT
data(:,3) = 2*data(:,3) -1;
nSubjects = size(data,1)/nTrials;
impossibleMat = nan(nSubjects, nImpTrials);
beforeImpossibleMat = nan(nSubjects, nImpTrials);
for nSub = 1:nSubjects
    matSubject = data( ( 1+(nSub-1)*nTrials ):(nSub*nTrials) , :);
    locImpossibleSubject =  find( matSubject(:,2) == 0);
    locPossibleSubject = locImpossibleSubject-1;
    impossibleMat(nSub,:) = matSubject( locImpossibleSubject, 3 )';
    beforeImpossibleMat(nSub,:) = matSubject( locPossibleSubject, 3 )';
end
pBiasVect = (sum(impossibleMat==1,2)/nImpTrials); % pUp in impossible trials 
pBeforeBiasVect = (sum(beforeImpossibleMat==1,2)/nImpTrials); % pUp in preceding trials 
biasVect = -1+2*(sum(impossibleMat==1,2)/nImpTrials);


% compute the PDF of pPrev using combinatorics: 
% (per-subjects impossible trials with previous trial)
forPprevImpPrev = (impossibleMat .* beforeImpossibleMat) == 1;
pPrevImpPrev_sub = sum(forPprevImpPrev,2)/ nImpTrials;
numUpImp_sub = sum( (impossibleMat == 1), 2);
numUpPos_sub = sum( (beforeImpossibleMat == 1), 2);
bySubjUpperLim = nan(nSubjects, 1);
bySubjLowerLim = nan(nSubjects, 1);
bySubjUpperLimPos = nan(nSubjects, 1);
bySubjLowerLimPos = nan(nSubjects, 1);
bySubj_pPrevPDF = nan( nSubjects, nImpTrials+1 );
load('prob4pPrevNotBiggerThanP2_data_SHORT.mat');
% note that significance is above or below the upper and lower lims,
% respectively.
for sub = 1:nSubjects
    bySubjUpperLim(sub) = ...
        CI90_okay_Up.(['impUp' num2str(numUpImp_sub(sub)) ...
        '_prevUp' num2str(numUpPos_sub(sub))]);
    bySubjLowerLim(sub) = ...
        CI90_okay_Down.(['impUp' num2str(numUpImp_sub(sub)) ...
        '_prevUp' num2str(numUpPos_sub(sub))]);
    bySubjUpperLimPos(sub) = ...
        CI_okay_UpPos.(['impUp' num2str(numUpImp_sub(sub)) ...
        '_prevUp' num2str(numUpPos_sub(sub))]);
    bySubjLowerLimPos(sub) = ...
        CI_okay_DownPos.(['impUp' num2str(numUpImp_sub(sub)) ...
        '_prevUp' num2str(numUpPos_sub(sub))]);
    bySubj_pPrevPDF(sub, : ) = PprevPDF.(['impUp' num2str(numUpImp_sub(sub)) ...
        '_prevUp' num2str(numUpPos_sub(sub))]);
end


% Plot Volume-like table:
sorted_pPrevImpPrev_sub = sortrows( [pPrevImpPrev_sub, bySubjLowerLim, ...
    bySubjLowerLimPos, bySubjUpperLim, bySubjUpperLimPos, pBiasVect, ...
    bySubj_pPrevPDF], 1:7 );
figure;
hold on;
plot( [sorted_pPrevImpPrev_sub(:,3), sorted_pPrevImpPrev_sub(:,5)]', ...
    [(1:nSubjects)' (1:nSubjects)']', 'g-', 'lineWidth', 1)
hold on;
plot( [-0.025+sorted_pPrevImpPrev_sub(:,2), +0.025+sorted_pPrevImpPrev_sub(:,4)]', ...
    [(1:nSubjects)' (1:nSubjects)']', 'k-', 'lineWidth', 1);
% we added 0.0125 for demonstration purposes, so that if one falls at the
% edge of the 90% condifence interval it will indeed show insignificance 
% on the volume-plot.
hold on;
plot( sorted_pPrevImpPrev_sub(:,1), 1:nSubjects, 'r.' ); hold on;
mean_pPrev = mean( sorted_pPrevImpPrev_sub(:,1) );
sem_pPrev = std( sorted_pPrevImpPrev_sub(:,1) ) / sqrt( -1 + length( sorted_pPrevImpPrev_sub(:,1) ) );
plot( [mean_pPrev, mean_pPrev], [-3 -1], 'r-', ...
    [mean_pPrev-sem_pPrev, mean_pPrev+sem_pPrev], [-2,-2], ...
    'r-', 'lineWidth', 1 );
xlim( [-0.05,1.05] ); ylim( [-0.075*nSubjects,1.02*nSubjects] );
xlabel('p_{repeat}'); ylabel('# participant')
xticks([0 0.25 0.5 0.75 1]); yticks( 0:20:nSubjects );
box on;
xlim([0,1]);


%% statistics for participant #100:

% pUp, ICB and pRepeat:
Participant100__pUp_ICB_pPrev = [sorted_pPrevImpPrev_sub(100,6), ...
    -1+2*sorted_pPrevImpPrev_sub(100,6), sorted_pPrevImpPrev_sub(100,1)]

% test if ICB is significantly different from chance:
pValue_pUp_participant100 = myBinomTest( ...
    nImpTrials*Participant100__pUp_ICB_pPrev(1), nImpTrials, 0.5 )

% compute the pValue correpsonding to the participant's pRepeat (using the 
% combinatorics-based PDF that was obtained by number of Up responses in
% the impossible and preceding possible trials):
pRepeatPDF_participant100 = sorted_pPrevImpPrev_sub(100,7:end);
loc_pPrevInPDF = find( (0:0.05:1) == Participant100__pUp_ICB_pPrev(3) );
pValue_pPrev_participant100 = sum( ...
    pRepeatPDF_participant100( loc_pPrevInPDF:end ) )

% If interested in verifying the obtained pValue, then go ahead and 
% uncomment the following:
%{
locParticipant = find( pPrevImpPrev_sub == max(pPrevImpPrev_sub) ); % participant no. in unsorted data
choices100Imp = impossibleMat(locParticipant,:); 
choices100Pos = beforeImpossibleMat(locParticipant,:);
nSim = 1e6;
pPrev_sim = nan(nSim,1);
for s = 1:nSim
    choices100Imp_sim = choices100Imp( randperm( nImpTrials ) ); % permute the order of impossible trials
    choices100Pos_sim = choices100Pos( randperm( nImpTrials ) ); % permute the order of vertical trials preceding the impossible
    forPprev_sim = ( (choices100Imp_sim .* choices100Pos_sim) == 1 ); 
    pPrev_sim(s) = sum(forPprev_sim,2)/ nImpTrials; % compute pRepeat for current permutation 
end
pValue_pPrev_participant100_sim = ...
    sum( pPrev_sim >= Participant100__pUp_ICB_pPrev(3) ) / nSim
%}
