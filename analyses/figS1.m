%% Supp Fig. 1

addpath('analyses','data');

data_BL_all = xlsread('BL_vertical_raw.xlsx', 1, 'C2:D12001');

nTrials = 120;
nSubjects = size(data_BL_all,1)/nTrials;

bias_imp_pos = zeros(nSubjects,3);
performance = nan(nSubjects,1);
allCon1Incon0Subs = [];
allCon1Incon0ICBs = [];

for i = 1:nSubjects
    sub_mat = data_BL_all( (1+nTrials*(i-1)):(nTrials*i), : );
    sub_mat_imp = sub_mat( sub_mat(:,1)==0, : ); % impossible trials
    polarBiasSub = ( mean(sub_mat_imp(:,2)) > 0.5 ) - ( mean(sub_mat_imp(:,2)) < 0.5 ) ;
    sub_mat_pos = sub_mat( sub_mat(:,1)~=0, : ); % possible trials
    sub_mat_ePos = sub_mat( ( ((2*sub_mat(:,2))-1) .* sub_mat(:,1) ) < 0, : ); % possible error trials
    if polarBiasSub~=0
        allCon1Incon0Subs = [allCon1Incon0Subs; ...
            polarBiasSub*(2*(sub_mat_ePos(:,2))-1)];
        allCon1Incon0ICBs = [allCon1Incon0ICBs; ...
            polarBiasSub*(2*(sub_mat_imp(:,2))-1)];
    end
    bias_imp_pos(i,1) = 2 * mean( sub_mat_imp(:,2) ) - 1;
    bias_imp_pos(i,2) = 2 * mean( sub_mat_pos(:,2) ) - 1; 
    bias_imp_pos(i,3) = 2 * mean( sub_mat_ePos(:,2) ) - 1; 
    performance(i) = mean(  ( ( -1 + 2 * (sub_mat_pos(:,1) > 0 ) ) .* ...
        ( -1 + 2 * (sub_mat_pos(:,2) > 0 ) ) )  >  0  );
end
bias_imp = bias_imp_pos(:,1);
bias_pos = bias_imp_pos(:,2);
bias_ePos = bias_imp_pos(:,3);
bias_imp_pos_unique = unique( [bias_imp, bias_pos], 'rows' );
all_mSize_check = nan(1, length(bias_imp_pos_unique) );

figure; 
% Impossible with possible trials:
for i = 1:length(bias_imp_pos_unique)
    mSize = sum( ( bias_imp_pos(:,1) == bias_imp_pos_unique(i,1) ) .* ...
        ( bias_imp_pos(:,2) == bias_imp_pos_unique(i,2) ) );
    all_mSize_check(i) = mSize;
    plotDots = plot( bias_imp_pos_unique(i,1), bias_imp_pos_unique(i,2), ...
        'Marker','o','MarkerSize',4+2*(mSize-1), ...
        'MarkerFaceColor', [0,0,0],'Color',[1 1 1]); hold on; 
end
% describe dot size:
nPairs = length(all_mSize_check)
nPairs_singleParticipants = sum(all_mSize_check==1)
nPairs_participantsSharingPairs = sum(all_mSize_check>1)
nParticipantsSharingPair = sum( all_mSize_check(all_mSize_check>1) )
unique_nParticipantsSharingPairs = unique(all_mSize_check);

% Orthogonal regression:
v = pca([bias_imp bias_pos]);
slope = v(2,1)/v(1,1);
k = mean( bias_pos ) - slope * mean( bias_imp );
plot( [-1,1], slope * [-1,1] + k, 'Color', [1,1,1], 'lineWidth', 2 ); hold on;
h = plot( [-1,1], slope * [-1,1] + k, 'Color', [0,0,0], 'lineWidth', 1 ); hold on;
mean_imp = mean( bias_imp );
mean_pos = mean( bias_pos );
sem_imp = std( bias_imp ) / sqrt( length(bias_imp) -1 );
sem_pos = std( bias_pos ) / sqrt( length(bias_pos) -1 );
[rho, pVal] = corr( bias_imp, bias_pos );
legend( [plotDots,h], ['data: rho = ' num2str(rho) ', p = ' num2str(pVal)], ...
    ['Ortho. reg: ICB_{pos} = ' num2str(slope) ' * ICB_{imp}'], ...
    'Location','SouthOutside');
xlabel('ICB impossible');
ylabel('ICB possible')
box off; axis square; ylim([-0.3,0.3]); xlim([-1,1]);
ggg = gca; ggg.XMinorTick = 'on'; ggg.YMinorTick = 'on';
xticks(-1:0.5:1); yticks(-0.3:0.15:0.3);    


%% Analysis of ICB in impossible and possible ERROR(!) trials:
v = pca([bias_imp bias_ePos]);
slope = v(2,1)/v(1,1);
[rho, pVal] = corr( bias_imp(~isnan(bias_ePos)), bias_ePos(~isnan(bias_ePos)) )
p_con_error = sum(allCon1Incon0Subs == 1) / length(allCon1Incon0Subs)
p_incon_error = sum(allCon1Incon0Subs == -1) / length(allCon1Incon0Subs)


%% Average performance +- std:
overAllPerformance = mean( performance )
stdPerformance = std( performance )
