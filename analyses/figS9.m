%% supp Fig. 9: The 'drift bias' DDM and 'IC+drift bias' DDM but not the 
% 'IC bias' DDM account well for decision-makers' responses. 


%% load PPC data. 

% For each dataset, we simulated responses (choices and reaction time) 
% using the corresponding posteriors obtained from the HDDM procedure 
% (2,000 for each participant and 10,000 responses for each recurrent network).

% Run PPC_loop function to run on datasets and models, unless output 
% data exists. This runs over observed and expected data (based on
% average posetior of each participant in the biased ‘IC+drift’ DDM .

addpath('analyses','data');

isPPC = exist('PPC_loop_data.mat'); %
if isPPC == 2
    load('PPC_loop_data.mat') 
else 
    PPC_loop
end

titleName = {'Vertical bisection', 'g = 2.7', 'g = 3', 'g = 3.4', 'Motor task'};
plotDataOrder = [1, 5, 4, 3, 2];
modelName_vect = {'none', 'drift', 'ic', 'both'}; 
dataName_vect = {'BL', 'g27', 'g3', 'g34', 'motor1', ...
    'motor2', 'motor3', 'motor4', 'motor5', 'motor6', 'motor7', ...
    'motor8', 'motor9', 'motor10'}; 
forModelsCols = [1 0 1; 0.0431 0.5020 0.2510; ...
    0.4980 0.1843 0.5451; 0 0 0; 0.93 0.69 0.13];
modelMarker = {'*', 'o', 's', 'x'};
markerSize = [4 4 4 6];
modelMarkerEdge = [0 0 0; 1 1 1; 1 1 1; 0 0 0];

figure;


%% supp Fig. 9a: Expected (simulated response) vs. observed probability of 
% choice for each observer:

for dataNum = 1:5
    
    tit = titleName{dataNum};
    dataName = dataName_vect(dataNum);
    subplot( 2, 5, plotDataOrder(dataNum) );
    
    % running over DDMs (bisaed IC, bisaed drifed, bisaed both):
    for modelNum = [3 2 4] 
        modelName = modelName_vect(modelNum);
        mark = modelMarker{modelNum};
        pE_pO_qE1to5_qO1to5 = [];
        
        % load relevant dataset-model data:
        if dataNum == 5 % if motor then merge data of all 10 pairs 
            for dataNumMotor = 5:14
                motorpE_pO_adding = ...
                    PPC_analysis_struct.( dataName_vect{dataNumMotor} ).( modelName_vect{modelNum} );
                pE_pO_qE1to5_qO1to5 = [pE_pO_qE1to5_qO1to5; motorpE_pO_adding]; 
            end
        else
            pE_pO_qE1to5_qO1to5 = ...
                PPC_analysis_struct.( dataName_vect{dataNum} ).( modelName_vect{modelNum} );
        end
        
        % plot observeda nd expected datapoints:
        plotEO = plot( pE_pO_qE1to5_qO1to5(:,2), pE_pO_qE1to5_qO1to5(:,1), ...
            'Marker', mark, 'LineStyle', 'none', 'markerSize', markerSize(modelNum), ...
            'MarkerFaceColor', forModelsCols(modelNum,:), ...
            'MarkerEdgeColor', modelMarkerEdge(modelNum,:) );
        
        % plot orthogonal regression (not shown in figures):
        v = pca([pE_pO_qE1to5_qO1to5(:,2) pE_pO_qE1to5_qO1to5(:,1)]); hold on;
        slope = v(2,1)/v(1,1);
        slope_vect_PPC(modelNum) = slope;
        k = mean( pE_pO_qE1to5_qO1to5(:,1) ) - slope * mean( pE_pO_qE1to5_qO1to5(:,2) );
        plot( [-1,1], slope * [-1,1] + k, 'Color', [1,1,1], 'lineWidth', 2 ); hold on;
        plotPPC(modelNum-1) = plot( [-1,1], slope * [-1,1] + k, ...
            'Color', forModelsCols(modelNum,:), 'lineWidth', 1 ); hold on;
    end
    
    xlim( [0,1] ); ylim( [0,1] );
    xlabel('Observed p'); ylabel('Expected p'); 
    box off; axis square;
    title( tit );
    
    % show plots of orthogonal regression:
    legend( [plotPPC(1), plotPPC(2), plotPPC(3)], ...
        {[modelName_vect{2} ', slope=' num2str(slope_vect_PPC(2))], ...
        [modelName_vect{3} ', slope=' num2str(slope_vect_PPC(3))], ...
        [modelName_vect{4} ', slope=' num2str(slope_vect_PPC(4))]}, ...
        'Location','SouthOutside' );
    
end 


%% FIG S9b: observed vs. expected (simulated response) RT distributions:

xLimData = [-5, 5; -5, 5; -1, 1; -0.5, 0.5; -5, 5];
yUpLim = [2, 4, 8, 12, 2];

for dataNum = [1 5 2 3 4]
    
    dataName = dataName_vect(dataNum);
    subplot( 2, 5, 5 + plotDataOrder(dataNum) );
    
    % running over DDMs (bisaed both, bisaed IC, bisaed drift):
    for modelNum = 4:-1:2
        RT_observed= [];
        RT_expected = [];
        
        % load relevant dataset-model data:
        if dataNum == 5 % if motor then merge data of all 10 pairs 
            for dataNumMotor = 5:14
                RT_observed_adding = PPC_RT_real.( dataName_vect{dataNumMotor} ).observed;
                RT_observed = [RT_observed; RT_observed_adding];
                RT_expected_adding = PPC_RT_real.( dataName_vect{dataNumMotor} ).( modelName_vect{modelNum} );
                RT_expected = [RT_expected; RT_expected_adding];
            end
        else
            RT_observed = PPC_RT_real.( dataName_vect{dataNum} ).observed;
            RT_expected = PPC_RT_real.( dataName_vect{dataNum} ).( modelName_vect{modelNum} );
        end
        
        % Determine whether plus/ minus ICB for positive/negative RTs:
        forPM_observed = (2 * sum( RT_observed > 0, 2 , 'omitnan') ./ ...
            sum( ~isnan(RT_observed), 2 ) ) - 1;
        % Disregard observers with ICB = 0:
        forPM_observed( forPM_observed == 0 ) = NaN;
        % Determine positive (decision congruent with ICB polarity) and 
        % negative (decision congruent with ICB polarity) RTs:
        RT_observed = ( 2 * ( forPM_observed > 0 ) - 1 ) .* RT_observed;
        RT_expected = ( 2 * ( forPM_observed > 0 ) - 1 ) .* RT_expected;
        RT_observed = reshape( RT_observed, [], 1 );
        RT_expected = reshape( RT_expected, [], 1 );
        % omit observers with ICB=0: 
        RT_observed = RT_observed( isnan(RT_observed) ~= 1 ); 
        RT_expected = RT_expected( isnan(RT_expected) ~= 1 );
        
        % Binning for positive and negative RTs:
        if (dataNum == 1) || (dataNum == 5)
            nBins = 32; 
            RT_observed = RT_observed(RT_observed<=3);
            RT_expected = RT_expected(RT_expected<=3);
        else
            nBins = 145;
        end
        
        binWidth = 0.025;
        max_edge = mean(abs(RT_observed)) + 7 * std(abs(RT_observed)); 
        
        % Plot congruent RTs:
        [hist_RT_observed, edges] = histcounts( RT_observed(RT_observed>=0), nBins, ...
            'Normalization', 'pdf', 'BinLimits',[0,max_edge]);
        hist_RT_expected = histcounts( RT_expected(RT_expected>=0), 'binEdges', edges, 'Normalization', 'pdf' );
        centers = 0.5*(edges(2)-edges(1)) + edges(1:(end-1)); 
        if modelNum == 4 % if 1st model then also plot observed
            a = stairs( centers, hist_RT_observed, ...
                'Color', [0.7 0.7 0.7], 'lineWidth', 2 ); hold on; 
        end
        xe(modelNum) = plot( centers, hist_RT_expected, ...
            'Color', forModelsCols(modelNum,:), 'lineWidth', 1 ); hold on;
        
        % Plot incongruent RTs:
        [hist_RT_observed, edges] = histcounts( RT_observed(RT_observed<0), nBins, ...
            'Normalization', 'pdf', 'BinLimits',[-max_edge,0]);
        hist_RT_expected = histcounts( RT_expected(RT_expected<0), 'binEdges', edges, 'Normalization', 'pdf' );
        centers = 0.5*(edges(2)-edges(1)) + edges(1:(end-1));
        if modelNum == 4 % if 1st model then also plot observed
            a = stairs( centers, hist_RT_observed, ...
                'Color', [0.7 0.7 0.7], 'lineWidth', 2 ); hold on; 
        end
        plot( centers, hist_RT_expected, ...
            'Color', forModelsCols(modelNum,:), 'lineWidth', 1 ); hold on;
        
    end
    
    xlabel('RT [s]'); ylabel('PDF'); axis square; box off;
    xlim( xLimData( dataNum, : ) ); ylim( [0, yUpLim( dataNum )] );
    
end
