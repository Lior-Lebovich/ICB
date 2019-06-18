%% supp Fig. 10: Plot fitted threshold/drift/IC vs. observed bias (BL,motor,I&F ).:


%% load PPC data - we will only use it for observed probabilities.
% Run PPC_loop function to run on datasets and models, unless output 
% data exists. 

addpath('analyses','data');

isPPC = exist('PPC_loop_data.mat'); 
if isPPC == 2
    load('PPC_loop_data.mat') 
else 
    PPC_loop
end

dataName_vect = {'BL', 'g27', 'g3', 'g34', ...
    'motor1', 'motor2', 'motor3', 'motor4', 'motor5', ...
    'motor6', 'motor7', 'motor8', 'motor9', 'motor10'}; 
titleName = {'Vertical bisection', 'g = 2.7', 'g = 3', 'g = 3.4', 'Motor task'};
PPC_link = 'hddm related code and data\';
threshLim = [1, 2.4; 1, 3; 0.7, 1.1; 0.55, 0.75; 1, 2.4];
driftLim = [0, 3; 0, 8; 0, 8; 0, 8; 0, 3];
plotDataOrder = [1, 5, 4, 3, 2];

figure;

% running over data (5):
for dataNum = 1:5
    dataName = dataName_vect{dataNum};
    stat_data_threshold = [];
    stat_data_drift = [];
    
    % #observers in dataset:
    if strcmp(dataName, 'BL')
        nSubjects = 100;
    elseif ( strcmp(dataName, 'g27') || strcmp(dataName, 'g3') || ...
        strcmp(dataName, 'g34') )
        nSubjects = 200;
    else
        nSubjects = 20;
    end
    
    % load observed p:
    pE_pO_qE1to5_qO1to5 = [];
    if dataNum == 5 % if motor then merge data of all 10 pairs 
        for dataNumMotor = 5:14
            motorpE_pO_adding = ...
            PPC_analysis_struct.( dataName_vect{dataNumMotor} ).drift;
            pE_pO_qE1to5_qO1to5 = [pE_pO_qE1to5_qO1to5; motorpE_pO_adding];
        end
        dataName2 = 'motor';
    else
        pE_pO_qE1to5_qO1to5 = ...
        PPC_analysis_struct.( dataName_vect{dataNum} ).drift;
        dataName2 = dataName;
    end
    p_observed = pE_pO_qE1to5_qO1to5(:,2);
    
    % Load threshold and drift posterior averages:
    if dataNum == 5 % if motor then merge data of all 10 pairs 
        for dataNumMotor = 5:14
            dataLink = PPC_link + string( dataName_vect{dataNumMotor} ) + ...
            '\informative priors 40K\drift\stats_csv_' + ...
            string(dataName2) + '_drift_no_tttv.csv';
            stat_data_threshold_adding =  csvread( dataLink ,3 ,1, [3, 1, 3+nSubjects-1, 1] );
            stat_data_threshold =[stat_data_threshold; stat_data_threshold_adding];
            stat_data_drift_adding =  csvread( dataLink, 2+3+nSubjects, 1, [2+3+nSubjects, 1, 2+3+2*nSubjects-1, 1] );
            stat_data_drift = [stat_data_drift; stat_data_drift_adding];
            for_stat_data_threshold_pop_adding =  csvread( dataLink ,1 ,1, [1, 1, 1, 1] );
            stat_data_threshold_pop_adding = for_stat_data_threshold_pop_adding * ones( nSubjects, 1 );
            stat_data_threshold_pop = [stat_data_threshold_pop; stat_data_threshold_pop_adding];
            for_stat_data_drift_pop_adding =  csvread( dataLink, 3+nSubjects, 1, [3+nSubjects, 1, 3+nSubjects, 1] );
            stat_data_drift_pop_adding = for_stat_data_drift_pop_adding * ones( nSubjects, 1 );
            stat_data_drift_pop = [stat_data_drift_pop; stat_data_drift_pop_adding];
            % for population average:
            stat_data_threshold_m_s_adding = csvread( dataLink ,1 ,1, [1, 1, 2, 1] );
            stat_data_threshold_m_s = ...
            [stat_data_threshold_m_s stat_data_threshold_m_s_adding];
            stat_data_drift_m_s_adding = ...
            csvread( dataLink, 3+nSubjects, 1, [3+nSubjects, 1, 4+nSubjects, 1] );
            stat_data_drift_m_s = ...
            [stat_data_drift_m_s stat_data_drift_m_s_adding];
        end
        stat_data_threshold_m_s = mean( stat_data_threshold_m_s, 2 );
        stat_data_drift_m_s = mean( stat_data_drift_m_s, 2 );
        else
        dataLink = PPC_link + string(dataName) + ...
        '\informative priors 40K\drift\stats_csv_' + ...
        string(dataName2) + '_drift_no_tttv.csv';
        stat_data_threshold =  csvread( dataLink ,3 ,1, [3, 1, 3+nSubjects-1, 1] );
        stat_data_drift =  csvread( dataLink, 2+3+nSubjects, 1, [2+3+nSubjects, 1, 2+3+2*nSubjects-1, 1] );
        for_stat_data_threshold_pop =  csvread( dataLink ,1 ,1, [1, 1, 1, 1] );
        stat_data_threshold_pop =  for_stat_data_threshold_pop * ones( nSubjects, 1 );
        for_stat_data_drift_pop =  csvread( dataLink, 3+nSubjects, 1, [3+nSubjects, 1, 3+nSubjects, 1] );
        stat_data_drift_pop = for_stat_data_drift_pop * ones( nSubjects, 1 );
        % for population average:
        stat_data_threshold_m_s = csvread( dataLink ,1 ,1, [1, 1, 2, 1] );
        stat_data_drift_m_s = csvread( dataLink, 3+nSubjects, 1, [3+nSubjects, 1, 4+nSubjects, 1] );
    end
    stat_data_drift_abs = abs( stat_data_drift ); % drift magnitude
    p_observed_max = max( [p_observed, 1-p_observed]' )';
    
    % plot observed p vs. average threshold:
    subplot(3,5,plotDataOrder(dataNum));
    plot( p_observed_max, stat_data_threshold, ...
    'Marker','h','LineStyle','none', 'MarkerSize', 4, ...
    'MarkerFaceColor', [.93 .69 .13],'Color',[1 1 1] );
    if dataNum == 1
        ylabel('Threshold (a)')
    end
    xlabel('Observed p_{bias}'); axis square;
    set(gca, 'FontName', 'Times'); box off;
    [rho_corr,p_corr] = corr( p_observed_max, stat_data_threshold, 'Type', 'Spearman'  );
    title( titleName{dataNum} );
    xlim( [0.5, 1] ); xticks( 0.5:0.25:1 );
    ylim( threshLim(dataNum,:) );
    
    % plot observed p vs. average drift:
    subplot(3,5,plotDataOrder(dataNum)+5);
    plot( p_observed_max, stat_data_drift_abs, 'Marker','o', ...
    'LineStyle','none', 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.0431 0.5020 0.2510],'Color',[1 1 1] );
    if dataNum == 1
        ylabel('Drift magnitude (|A|)')
    end
    xlabel('Observed p_{bias}'); xlim([.5,1]); xticks([0.5,0.75,1]); axis square;
    set(gca, 'FontName', 'Times'); box off;
    xlim( [0.5, 1] ); xticks( 0.5:0.25:1 );
    ylim( driftLim(dataNum,:) );
    
    % Compute p_expected 
    subplot(3,5,plotDataOrder(dataNum)+10);
    p_given_mjDrift_mjThresh = 1 ./ ... % observer mean drift and threshold
        ( 1 + exp( -1 * stat_data_drift_abs .* stat_data_threshold ) );
    p_given_mjDrift_mPopThresh = 1 ./ ... % observer mean drift, pop. average threshold
        ( 1 + exp( -1 * stat_data_drift_abs .* mean(stat_data_threshold) ) );
    p_given_mjThresh_mPopDrift = 1 ./ ... % observer mean threshold, pop. average drift
        ( 1 + exp( -1 * mean(stat_data_drift_abs) .* stat_data_threshold ) );
    
    % plot pO vs. pE given observer drift and pop threshold:
    pSub_drift = plot( p_observed_max, p_given_mjDrift_mPopThresh, ...
        'Marker','o','LineStyle','none', 'MarkerFaceColor',...
        [0.0431 0.5020 0.2510],'Color',[1 1 1], 'MarkerSize', 4 );
    hold on;
    % orthogonal regression:
    v = pca([p_observed_max p_given_mjDrift_mPopThresh]);
    slope_mjDrift_mPopThresh = v(2,1)/v(1,1);
    k = mean( p_given_mjDrift_mPopThresh ) - slope_mjDrift_mPopThresh * mean( p_observed_max );
    plot( [0.5,1], slope_mjDrift_mPopThresh * [0.5,1] + k, 'Color', [1,1,1], 'lineWidth', 2 ); hold on;
    pSub_drift_reg = plot( [0.5,1], slope_mjDrift_mPopThresh * [0.5,1] + k, 'Color', [0.0431 0.5020 0.2510], 'lineWidth', 1 ); hold on;
    
    % plot pO vs. pE given observer threshold and pop drift:
    pSub_thresh = plot( p_observed_max, p_given_mjThresh_mPopDrift, ...
        'Marker','h','LineStyle','none', 'MarkerFaceColor',[.93 .69 .13],...
        'MarkerSize', 4, 'Color',[1 1 1] ); 
    hold on;
    % orthogonal regression:
    v = pca([p_observed_max p_given_mjThresh_mPopDrift]);
    slope_mjThresh_mPopDrift = v(2,1)/v(1,1);
    k = mean( p_given_mjThresh_mPopDrift ) - slope_mjThresh_mPopDrift * mean( p_observed_max );
    plot( [0.5,1], slope_mjThresh_mPopDrift * [0.5,1] + k, 'Color', [1,1,1], 'lineWidth', 2 ); hold on;
    pSub_thresh_reg = plot( [0.5,1], slope_mjThresh_mPopDrift * [0.5,1] + k, 'Color', [.93 .69 .13], 'lineWidth', 1 ); hold on;
    
    % plot pO vs. pE given observer drift and threshold:
    pSub = plot( p_observed_max, p_given_mjDrift_mjThresh, ...
        'k', 'lineStyle', 'none', 'marker', 'x', 'MarkerEdgeColor', 'k' );
    
    xlabel('Observed p_{bias}');
    xlim([.5,1]); ylim([.5,1]);
    xticks([0.5,0.75,1]); yticks([0.5,0.75,1]); axis square; box off;
    set(gca, 'FontName', 'Times');
    if dataNum == 1
        ylabel('Expected p_{bias}');
    end
    legend( [pSub_thresh_reg, pSub_drift_reg], ...
    {['slope=' num2str( 0.01 * round( 100 * slope_mjThresh_mPopDrift ) )], ...
    ['slope=' num2str( 0.01 * round( 100 * slope_mjDrift_mPopThresh ) )]}, ...
    'Location','SouthEast' );

end
