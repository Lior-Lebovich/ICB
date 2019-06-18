%% Fig. 1d, Fig. 2d, supp Fig. 7: Relative contributions of the drift-bias 
% and IC-bias to the ICBs in the 'IC+drift bias' DDM. 
% Computing observed vs. expected pBias, based on average posteriors of 
% each participant in the biased ‘IC+drift’ DDM.


%% load PPC data - we will only use it for observed probabilities.
% Run PPC_loop function to run on datasets and models, unless output 
% data exists. 

addpath('analyses','data');

isPPC = exist('PPC_loop_data.mat'); %
if isPPC == 2
    load('PPC_loop_data.mat') 
else 
    PPC_loop
end

dataName_vect = {'BL', 'g27', 'g3', 'g34', ...
    'motor1', 'motor2', 'motor3', 'motor4', 'motor5', ...
    'motor6', 'motor7', 'motor8', 'motor9', 'motor10'}; 
titleName = {'Vertical bisection', 'g = 2.7', 'g = 3', 'g = 3.4', 'Motor task'};
modelName_vect = {'none', 'drift', 'ic', 'both'}; 
forModelsCols = [1 0 1; 0.0431 0.5020 0.2510; ...
    0.4980 0.1843 0.5451; 0 0 0; 0.93 0.69 0.13];
modelsNameForDIC = {'b. drift','c. IC','a. drift + IC'};
modelMarker = {'*', 'o', 's', 'x'};
markerSize = [4 4 4 6];
modelMarkerEdge = [0 0 0; 1 1 1; 1 1 1; 0 0 0];
nModelsInclude = 4;
plotDataOrder = [1, 5, 4, 3, 2];

figure;


%% Relative contributions of the drift-bias 
 
for dataNum = 1:5
    
    dataName = dataName_vect{dataNum};
    tit = titleName{dataNum};
    subplot( 1, 5, plotDataOrder(dataNum) );
    
    % # observers in a dataset:
    if strcmp(dataName, 'BL') 
        nSubjects = 100;
    elseif ( strcmp(dataName, 'g27') || strcmp(dataName, 'g3') || ...
            strcmp(dataName, 'g34') )
        nSubjects = 200;
    else
        nSubjects = 20;
    end
    
    % compute p observed:
    pE_pO_qE1to5_qO1to5 = [];
    if dataNum == 5 % if motor then merge data of all 10 pairs 
        for dataNumMotor = 5:14
            motorpE_pO_adding = ...
                PPC_analysis_struct.( dataName_vect{dataNumMotor} ).both;
            pE_pO_qE1to5_qO1to5 = [pE_pO_qE1to5_qO1to5; motorpE_pO_adding]; 
        end
        dataName2 = 'motor';
    else
        pE_pO_qE1to5_qO1to5 = ...
            PPC_analysis_struct.( dataName_vect{dataNum} ).both;
        dataName2 = dataName;
    end
    p_observed = pE_pO_qE1to5_qO1to5(:,2);
    
    % Get poterior averaged drift, IC and threshold:
    stat_data_threshold = [];
    stat_data_drift = [];
    stat_data_ic = [];
    if dataNum == 5 % if motor then merge data of all 10 pairs 
        for dataNumMotor = 5:14 
            dataPath = 'hddm related code and data\' + string( dataName_vect{dataNumMotor} ) + ...
                '\informative priors 40K\both\stats_csv_' + ...
                string(dataName2) + '_both_no_tttv.csv';
            stat_data_threshold_adding =  csvread( dataPath, ...
                1+(1*2)+(nSubjects*(1-1)), 1, ...
                [1+(1*2)+(nSubjects*(1-1)), 1, (2*1)+(nSubjects*1), 1] );
            stat_data_threshold = [stat_data_threshold; stat_data_threshold_adding];
            stat_data_drift_adding =  csvread( dataPath, ...
                1+(2*2)+(nSubjects*(2-1)), 1, ...
                [1+(2*2)+(nSubjects*(2-1)), 1, (2*2)+(nSubjects*2), 1] );
            stat_data_drift = [stat_data_drift; stat_data_drift_adding];
            stat_data_ic_adding =  csvread( dataPath, ...
                1+(2*4)+(nSubjects*(4-1)), 1, ...
                [1+(2*4)+(nSubjects*(4-1)), 1, (2*4)+(nSubjects*4), 1] );
            stat_data_ic =[stat_data_ic; stat_data_ic_adding];
        end
    else
        dataPath = 'hddm related code and data\' + string(dataName) + ...
            '\informative priors 40K\both\stats_csv_' + ...
            string(dataName2) + '_both_no_tttv.csv';
        stat_data_threshold =  csvread( dataPath, ...
            1+(1*2)+(nSubjects*(1-1)), 1, ...
            [1+(1*2)+(nSubjects*(1-1)), 1, (2*1)+(nSubjects*1), 1] );
        stat_data_drift =  csvread( dataPath, ...
            1+(2*2)+(nSubjects*(2-1)), 1, ...
            [1+(2*2)+(nSubjects*(2-1)), 1, (2*2)+(nSubjects*2), 1] );
        stat_data_ic =  csvread( dataPath , ...
            1+(2*4)+(nSubjects*(4-1)), 1, ...
            [1+(2*4)+(nSubjects*(4-1)), 1, (2*4)+(nSubjects*4), 1] );
    end
    
    % Compute expected p, based on (by-observer) average postriors:
    p_given_mean_paramters_drift = ...
        ( 1 + exp( -1 * stat_data_drift .* stat_data_threshold ) ) .^ -1;
    p_given_mean_paramters_ic = stat_data_ic;
    p_given_mean_paramters_all = ...
        ( ( 1 + exp( -1 * stat_data_drift .* stat_data_threshold ) ) .^ -1 ) + ...
        ( ( 1 - exp( -1 * stat_data_drift .* stat_data_threshold .* ( -1 + 2 * stat_data_ic ) ) ) ./ ...
        ( exp( +1 * stat_data_drift .* stat_data_threshold ) - ...
        exp( -1 * stat_data_drift .* stat_data_threshold ) ) );
    
    % plot p observed vs. expected given unbiased drift and biased IC:
    plot( p_observed, p_given_mean_paramters_ic, ...
        'Marker', modelMarker{3}, 'LineStyle' ,'none', 'MarkerFaceColor',...
        forModelsCols(3,:), 'markerSize', markerSize(3), ...
        'MarkerEdgeColor', modelMarkerEdge(3,:) ); 
    hold on;
    % plot orthoginal regression:
    v = pca([p_observed p_given_mean_paramters_ic]);
    beta_ic = v(2,1)/v(1,1);
    k = mean( p_given_mean_paramters_ic ) - beta_ic * mean( p_observed );
    plot( [0,1], beta_ic * [0,1] + k, 'Color', [1,1,1], 'lineWidth', 2 ); hold on;
    h_ic = plot( [0,1], beta_ic * [0,1] + k, 'Color', forModelsCols(3,:), 'lineWidth', 1 ); hold on;
    
    % plot p observed vs. expected given biased drift and unbiased IC:
    plot( p_observed, p_given_mean_paramters_drift, ...
        'Marker', modelMarker{2}, 'LineStyle','none', 'MarkerFaceColor',...
        forModelsCols(2,:), 'markerSize', markerSize(2), ...
        'MarkerEdgeColor', modelMarkerEdge(2,:) ); 
    hold on;
    % plot orthoginal regression:
    v = pca([p_observed p_given_mean_paramters_drift]);
    beta_drift = v(2,1)/v(1,1);
    k = mean( p_given_mean_paramters_drift ) - beta_drift * mean( p_observed );
    plot( [0,1], beta_drift * [0,1] + k, 'Color', [1,1,1], 'lineWidth', 2 ); hold on;
    h_drift = plot( [0,1], beta_drift * [0,1] + k, 'Color', forModelsCols(2,:), 'lineWidth', 1 ); hold on;
    
    % plot p observed vs. expected given biased drift and biased IC:
    plot( p_observed, p_given_mean_paramters_all, ...
        'Marker',modelMarker{4},'LineStyle','none', 'MarkerFaceColor',...
        forModelsCols(4,:), 'markerSize', markerSize(4), ...
        'MarkerEdgeColor', modelMarkerEdge(4,:) ); 
    hold on;
    % plot orthoginal regression:
    v = pca([p_observed p_given_mean_paramters_all]);
    beta_both = v(2,1)/v(1,1);
    k = mean( p_given_mean_paramters_all ) - beta_both * mean( p_observed );
    plot( [0,1], beta_both * [0,1] + k, 'Color', [1,1,1], 'lineWidth', 2 ); hold on;
    h_both = plot( [0,1], beta_both * [0,1] + k, 'Color', forModelsCols(4,:), 'lineWidth', 1 ); hold on;
    
    if dataNum == 1
        ylabel('Expected p')
    end
    xlabel('Observed p');
    title( tit );
    xlim( [0,1] ); ylim( [0,1] );
    legend( [h_both, h_drift, h_ic], ...
        { ['slope = ' num2str(beta_both)], ...
        ['slope = ' num2str(beta_drift)], ...
        ['slope = ' num2str(beta_ic)] }, ...
        'Location','SouthOutside' );
    axis square; box off;
    xticks(0:0.5:1); yticks(0:0.5:1);

end
