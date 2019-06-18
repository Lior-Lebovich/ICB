%% supp Fig. S8:

addpath('analyses','data');

% data in the 4 files below correspond to pUp of 200 networks, each 
% simulated over 500 trials.
fileName = {'recurrent_g3_original', 'recurrent_g3_nE20K_nI20K.mat', ...
    'recurrent_g3_2N.mat', 'recurrent_g3_factor2interactions.mat' }; 
dataName = {'Original', 'nE = nI = 20,000', '2N', 'interactions' };

nbins = 25;
edges500_edges = linspace(0,1,26); %%(1/500)*[-0.001:20:239.999,260.001:20:500.001];
edges500_preference = (-1+(1/nbins)):(2/nbins):(1-(1/nbins));

% 10a. original data;   10b. 2N;   10c. ratio nE = nI = 20,000;  10d.
% factor 2 interactions.
figure; 

for j = 1:4
    
    % load data:
    pUp = load( fileName{j} );
    p = pUp.data;
    stdv = std(-1+2*p);
        
    % plot:
    subplot(2,2,j);
    [h, edges] = histcounts( p, 'Normalization', 'pdf', 'binEdges', edges500_edges );
    %histogram( p, 'Normalization', 'pdf', 'binEdges', edges500_edges, ...
    %    'FaceColor',[.5,.5,.5],'edgeColor',[1 1 1] );
    bar( edges500_preference, 0.5*h, 'FaceColor',[.5,.5,.5],'edgeColor','none' );
    xlabel('ICB'); ylabel('PDF'); box off; ylim([0,1.5]); xlim([-1,1]);
    title( [dataName{j}, ', std = ' num2str(stdv)]  ); 
end
