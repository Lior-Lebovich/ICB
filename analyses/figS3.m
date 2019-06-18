%% supp. Fig. S3: distributions of choice biases for all pairs in the 
% motor task:

addpath('analyses','data');
load('motor_raw');
figure;

pair_bars_count = nan(10,21);
i = 0;
for pair = [7 6 2 10 8 9 5 4 1 3]
    i = i + 1;
    % compute hist counts for the pair:
    pair_bars_count(pair,:) = histcounts( -1 + 2 * mean( choices_all(:,pair,:) ), ...
        'binEdges',linspace(-1,1,22),'Normalization','count');
    % plot:
    subplot(5,2,i);
    bar( -1:0.1:1, pair_bars_count(pair,:), ...
        'FaceColor', (1/255) * colors(pair,:), ...
        'EdgeColor', 'none' );
    xlim([-1.05,1.05]); ylim([0,7]); box off;
    yticks(0:2:6);
    ylabel('# participants');
    xlabel('ICB'); 
    xticks(-1:1:1);
end