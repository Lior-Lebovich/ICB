%% supp Fig. 6: DDM analysis of the strength of selective inhibition in the
% recurrent network model - contribution of the threshold and drift.

addpath('analyses'); addpath('data');


%% Loading per network posterior-averaged threshold and drift:

nNetworks = 200;
stat_data_threshold = nan(nNetworks,3);
stat_data_drift = nan(nNetworks,3);
gNames = {'27', '3', '34'};

for dataNum = 1:3 % running over data (g=2.7, g=3, g=3.4)
    dataName = gNames{dataNum};
    dataPath = 'hddm related code and data\g' + string(dataName) + ...
        '\informative priors 40K\drift\stats_csv_g' + ...
        string(dataName) + '_drift_no_tttv.csv';
    stat_data_threshold(:,dataNum) =  csvread( dataPath ,3 ,1, [3, 1, 3+nNetworks-1, 1] );
    stat_data_drift(:,dataNum) =  csvread( dataPath, 2+3+nNetworks, 1, [2+3+nNetworks, 1, 2+3+2*nNetworks-1, 1] );
end


%% plot per network posterior-averaged threshold and drift:

stat_data_drift_absNew = abs( stat_data_drift ); 
for net = 1:nNetworks
    col = rand(1,3);
    subplot(1,2,1); % threshold
    plot( [1,2,3], stat_data_threshold(net,:), 'Color', col, 'lineWidth', 1 ); hold on;
    subplot(1,2,2); % drift
    plot( [1,2,3], stat_data_drift_absNew(net,:), 'Color', col, 'lineWidth', 1 ); hold on;
end


%% plot posterior-averaged threshold and drift, averaged over networks:

% threshold:
subplot(1,2,1); 
plot( 1:3, mean(stat_data_threshold), 'Color', [0, 0, 0], 'lineWidth', 4 );
errorbar( 1:3, mean(stat_data_threshold), ...
    std(stat_data_threshold) / sqrt( size(stat_data_threshold,1) ), ...
    'lineStyle', 'none', 'lineWidth', 2, 'Color', 'k' );
plot( 1:3, mean(stat_data_threshold), 'lineStyle', 'none', ....
    'marker', 'o', 'MarkerSize', 4, 'lineWidth', 0.5, ...
    'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', [0,0,0], 'Color', 'k' );
xticks(1:3); xticklabels({'g = 2.7','g = 3.0','g = 3.4'}); 
ylabel('Threshold (a)'); box off; axis square;

% drift:
subplot(1,2,2); 
plot( 1:3, mean(stat_data_drift_absNew), 'Color', [0, 0, 0], 'lineWidth', 4 );
errorbar( 1:3, mean(stat_data_drift_absNew), ...
    std(stat_data_drift_absNew) / sqrt( size(stat_data_drift_absNew,1) ), ...
    'lineStyle', 'none', 'lineWidth', 2, 'Color', 'k' );
plot( 1:3, mean(stat_data_drift_absNew), 'lineStyle', 'none', ...
    'marker', 'o', 'MarkerSize', 4, 'lineWidth', 0.5, ...
    'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', [0,0,0], 'Color', 'k' );
xticks(1:3); xticklabels({'g = 2.7','g = 3.0','g = 3.4'}); 
ylabel('Drift magnitude (|A|)'); box off; axis square;


%% Compute the average change from g=2.7 to g=3.4 for each averaged 
% posterior, averaged over networks:

avgChangeThresh_g27_g34 = -1 + mean( stat_data_threshold(:,3) ) ./ ...
    mean( stat_data_threshold(:,1) )
avgChangeDrift_g27_g34 = -1 + mean( stat_data_drift_absNew(:,3) ) ./ ...
    mean( stat_data_drift_absNew(:,1) )
