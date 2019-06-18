%% Figs. 1a + 1b:

addpath('analyses','data');


%% Fig. 1a: psychometric curves of 3 participants in the vertical bisection task

% read data from xlsx file of vertical bisection task: 
data_dev_isUp = xlsread('BL_vertical_raw.xlsx', 1, 'C2:D12001');

% compute pUp for each partcipant (row) and each deviation (col):
nTrials = 120;
nSubjects = size(data_dev_isUp,1)/nTrials;
dev = unique( data_dev_isUp(:,1) );
pUpByDl = nan( nSubjects, length(dev) );
nTrialsByDl = nan( nSubjects, length(dev) );
for i = 1:nSubjects
    dataSubj = data_dev_isUp( (1+120*(i-1)):(120*i), : );
    for j = 1:length(dev)
        choicesSubjDl = dataSubj( dataSubj(:,1)==dev(j), 2 );
        pUpByDl(i,j) = mean(choicesSubjDl);
        nTrialsByDl(i,j) = length(choicesSubjDl);
    end
end

% fit psychometric functions:        
ft = fittype( '(1+exp(-a*(x-b)))^(-1)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.7791 0.8427];
dlOverL = (1/200) * 2 * dev'; % relative offse = DeltaL/L
st_err = sqrt( ( pUpByDl .* (1-pUpByDl) ) ./ nTrialsByDl ); % SEM of pUp
subNumPsy = [33, 40, 79]; % representative participants
cols = {'blue', 'black', 'red'};
marks = {'o', 's', 'd'};

% plot:
figure;
for k = 1:3
    col = cols{k};
    mark = marks{k};
    f1B_curve = pUpByDl(subNumPsy(k),:);
    [xDataSS, yDataSS1] = prepareCurveData(dlOverL,f1B_curve);
    [fitresult, ~] = fit( xDataSS, yDataSS1, ft, opts );
    errorbar(dlOverL,pUpByDl(subNumPsy(k),:),st_err(subNumPsy(k),:),...
        'MarkerEdgeColor','none','MarkerSize',5,'Marker',mark,...
        'LineStyle','none','lineWidth',0.51,'color',col,'MarkerFaceColor',col); 
    hold on; 
    plot(fitresult,col); hold on;
end
xlim([min(dlOverL),max(dlOverL)]); ylim([0,1]); box off; legend off;
xlabel('\DeltaL/L'); ylabel('p_{Up}');
ggg = gca; ggg.XMinorTick = 'on'; ggg.YMinorTick = 'on';
xticks(-0.1:0.1:0.1); yticks(0:0.5:1);



%% Fig. 1b: ICB distribution for impossible trials of the vertical bisection task

locImpossible = find( dev == 0 ); % impossible trials
ICB_BL = -1 + 2 * pUpByDl( :, locImpossible ); % ICB = pUp - pDown 
edges = linspace(-1,1,22);
ICB_BL_pdf = histcounts( ICB_BL, 'Normalization', 'pdf', 'binEdges', edges );

% plot PDF:
figure;
bar( -1:0.1:1, ICB_BL_pdf, 'FaceColor', [.5 .5 .5], 'edgeColor', 'none' ); 
xlim([-1.05,1.05]); box off; 
xlabel('ICB'); ylabel('PDF');
ggg = gca; ggg.XMinorTick = 'on'; ggg.YMinorTick = 'on';
xticks(-1:1:1); yticks(0:0.5:1); hold on;

% plot ICBs that correspond to the 3 psychometric curves: 
for k = 1:3
    col = cols{k};
    PuP_subPsy = pUpByDl( subNumPsy(k), locImpossible );
    ICB_subPsy = -1 + 2 * pUpByDl( subNumPsy(k), locImpossible );
    relatedPDF = ICB_BL_pdf( find( 0:.05:1 == PuP_subPsy ) );
    text( ICB_subPsy, relatedPDF, '\downarrow', 'color', col, 'FontSize',15, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold');
    hold on;
end


%% Tests:

% Significant biases: Binomial tests, not corrected for multiple
% comparispons.
% Compute #significant, 2-sided binomial test:
sumSigDownBias = ...
    sum( myBinomTest( 20*pUpByDl( :, locImpossible ), 20, 1/2) < 0.05 );
% significant biases correspond to pUp<=0.25 or pUp>=0.75 (0.5<=|ICB|<=1):
pdfBinom =  pdf('Binomial',0:20,20,0.5);
maximalSigAlpha = sum(pdfBinom(1:6))*2; 
% Compute #significant Up/Down ICBs, 2-sided binomial test:
sumSigDownBias =  sum(ICB_BL <= -0.5); % # significant 'Down' ICBs
sumSigUpBias =  sum(ICB_BL >= 0.5); % # significant 'Up' ICBs
% Compute significance levels for the 3 participants in Fig. 1a:
pUp_sig_mat = nan(3,2);
for k =1:3
    pUp_sig_mat(k,1) = pUpByDl( subNumPsy(k), locImpossible );
    pUp_sig_mat(k,2) = myBinomTest( 20 * pUp_sig_mat(k,1), 20, 1/2 );
end

% bootstrap test for global bias:
nSim = 1e6;
nImpossibleTrials = 20;
sim_avgPup = nan( nSim,1 );
for sim = 1:nSim
    pUp_sim = (1 / nImpossibleTrials) * binornd( nImpossibleTrials, ...
        datasample( pUpByDl( :, locImpossible ), nSubjects )  );
    sim_avgPup(sim) = mean( pUp_sim );
end
avgPup_95CI = quantile( sim_avgPup, [.025, 0.975] )
    
% bootstrap the std:
%real_avgPup = 0.5; % fair Bernoulli process
real_avgPup = mean( pUpByDl( :, locImpossible ) ); %Bernoulli process with p = average pUp
real_stdPup = std( pUpByDl( :, locImpossible ) );
sim_stdPup = std( (1 / nImpossibleTrials) * ...
    binornd( nImpossibleTrials, real_avgPup, nSubjects, nSim  ) );
sigLevel = sum( sim_stdPup > real_stdPup ) / nSim;
