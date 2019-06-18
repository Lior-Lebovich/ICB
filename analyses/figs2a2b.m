%% Fig. 2a, Fig. 2b: ICBs in the motor task

addpath('analyses','data');
load('motor_raw');
figure;


%% Fig. 2a: distribution of ICBs of all participants (n= 20), for exmaple 
% pair of dots (pair 7):

% compute hist counts:
example_pair = 7;
bars_count_example = histcounts( -1 + 2 * mean( choices_all(:,example_pair,:) ), ...
    'binEdges',linspace(-1,1,22),'Normalization','count');
% plot:
subplot(2,1,1);
bar( -1:0.1:1, bars_count_example, ...
    'FaceColor', (1/255) * colors(example_pair,:), ...
    'EdgeColor', 'none' );
xlim([-1.05,1.05]); ylim([0,4]); box off;
xticks(-1:1:1); yticks(0:2:6);
xlabel('ICB'); ylabel('# participants');


%% Fig. 2b: distribution of ICBs for all 10 pairs of dots in the experiment

% compute the PDF for all participant x pais (20x10):
pair_bars_PDF = nan(10,21);
newCols = nan(10,3);
i = 0;
for pair = [1 2 3 5 4 7 8 6 9 10] 
    i = i + 1;
    pair_bars_PDF(i,:) = histcounts( -1 + 2 * mean( choices_all(:,pair,:) ), ...
        'binEdges',linspace(-1,1,22),'Normalization','pdf');
    newCols(i,:) = (1/255) * colors(pair,:);
end

% plot:
subplot(2,1,2);
b = bar( -1:0.1:1, 0.1 * (pair_bars_PDF'), 'Stacked' );
for pair = 1:10
    b(pair).FaceColor = newCols(pair,:);
    b(pair).EdgeColor = 'none';
end
xlabel('ICB'); ylabel('PDF');
xlim([-1.05,1.05]); ylim([0,1.1]); box off;
