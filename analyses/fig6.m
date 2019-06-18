%% Fig. 6: ICBs in the recurrent spiking network model


addpath('analyses'); addpath('data');

figure;
load('dataFig6.mat');
 
nTrials = 500;
nNetworks = 200;

psycho_edge = 0.135;

% For ICB PDF: because there are 501 possible p, we define the binEdes so 
% that there are 21 p's in the central bin and 20 p's in all other bins. 
p500_edges = (1/500)*[-0.001:20:239.999,260.001:20:500.001]; 
nbins = length( p500_edges ) - 1;
ICB500_centers = (-1+(1/nbins)):(2/nbins):(1-(1/nbins));

gNames = {'34', '30', '27'};
netNames = {'net05', 'net14', 'net17'};
netColors = {'black', 'blue', 'red'};
netMarkers = {'o', 's', 'd'};
gLoop = 0;


for g = [3.4, 3, 2.7] % columns of Fig. 6
    gLoop = gLoop + 1;
    netLoop = 0;
    gName = gNames{gLoop};
    
    
    % Fig. 6b: Distribution of ICBs for each g (200 networks):
    
    subplot( 3, 3, 3 + gLoop );
    
    % compute std(ICB) for bootstrap test (below):
    std_ICB.( ['g' gName] ) = std( -1 + 2 * pUpImpDat.( ['g' gName] ) );
    
    % compute the PDF, from the vector of pUp - 200 network, 500 trial each
    pUp_PDF = histcounts( pUpImpDat.( ['g' gName] ), ...
        'binEdges', p500_edges, 'Normalization', 'pdf' );
    ICB_PDF = 0.5 * pUp_PDF; 
    
    % plot:
    bar(ICB500_centers, ICB_PDF, 'FaceColor', [.5 .5 .5], 'edgeColor', 'none'); 
    xlim( [-1, 1] ); ylim( [0, 2.5] ); box off; 
    ff = gca; ff.YMinorTick = 'on';
    xticks( -1:1:1 ); yticks( [0, 2.5] );
    
    
    % Fig. 6a: psychometric curves for each g:
    
    for net = [5 ,14 ,17] % psychometric curves for each network in g
        netLoop  = netLoop +1;
        netColor = netColors{netLoop};
        netName = netNames{netLoop};
        dl_p = psychoDat.( ['g' gName] ).( netName );
        
        % plot psychometric curves:
        subplot( 3, 3, gLoop );
        errorbar(dl_p(:,1), dl_p(:,2),...
            sqrt( ( dl_p(:,2) .* (1-dl_p(:,2) ) ) / nTrials ),...
            'MarkerEdgeColor','none',...
            'MarkerSize',5,'Marker','square','LineStyle','-','lineWidth',1,...
            'color',netColor,'MarkerFaceColor',netColor,'MarkerSize',4);
        hold on;
        
        % for arrows in Fig. 6b, compute pUp in impossible trials:
        Psycho_pUpImp = dl_p( dl_p(:,1) == 0, 2 );
        loc_ICB = max( find( p500_edges <= Psycho_pUpImp) ) ;
        x_Psycho_ICB = ICB500_centers( loc_ICB );
        y_Psycho_ICB = ICB_PDF( loc_ICB );
        % plot network ICB on ICB PDF:
        subplot( 3, 3, 3 + gLoop );
        text( x_Psycho_ICB, y_Psycho_ICB, '\downarrow', ...
            'color', netColor, 'FontSize',15, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        hold on;
    end
    subplot( 3, 3, gLoop ); 
    box off; 
    xlim([-psycho_edge*1.05,psycho_edge*1.05]); ylim([-0.05,1]);
    title(['g = ' num2str(g)]);
    ff = gca; ff.YMinorTick = 'on';
    xticks( -.1:.1:.1 ); yticks( 0:0.5:1 );
    
    
    % Fig. 6c: Distribution of reaction times for each g (200 networks):
    
    subplot( 3, 3, 6 + gLoop );
    
    % compute the PDF of decision times:
    widthDT = binWidthDT.( ['g' gName] );
    [DT_PDF_y, edgesDT] = histcounts( DecTimeImpDat.( ['g' gName] ), ...
        'binWidth', widthDT, 'BinLimits', [0,6], 'Normalization' ,'pdf' );
    DT_PDF_x = edgesDT( 1:(end-1) ) + 0.5 * widthDT;
    % start PDF from DT=zero (otherwise it is the center of 1st bin):
    add_left = 1:widthDT:( DT_PDF_x(1) - widthDT );
    % end PDF from DT=3 (otherwise it is the center of last bin):
    add_right = ( DT_PDF_x(end) + widthDT ):widthDT:3.1;

    % plot DT PDF:
    area( [add_left, DT_PDF_x, add_right],...
        [zeros(size(add_left)), DT_PDF_y, zeros(size(add_right))], ...
        'FaceColor', [.5 .5 .5], 'EdgeColor', 'none' );
    box off; 
    xlabel('RT [s]'); ylabel('PDF');
    xlim( [0, 3] ); ylim( [0, yUpLimDT.( ['g' gName] )] );

end


%% Bootstrap test for the std of ICB distributions in fig 6b:

% bootsrapping a fair bernoulli process and computing the std(ICB):
nSim = 1e6;
bootFairBernoulli_stdICB = zeros(1,nSim);
for kk = 1:nSim
    bootFairBernoulli_stdICB(kk) = ...
        std( 2 * mean( round( rand( nTrials, nNetworks ) ) ) - 1 );
end

% compute the p.Value of each g:
gLoop = 0;
for g = [3.4, 3, 2.7]
    gLoop = gLoop + 1;
    gName = gNames{gLoop};
    pValsigICB.( ['g' gName] ) = ...
        sum( bootFairBernoulli_stdICB >= std_ICB.( ['g' gName] ) ) / nSim;
end
pValsigICB
