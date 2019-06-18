%% Fig. 7: Conditional bias functions (CBF)

addpath('analyses','data');

% Fig. 7a: CBF for Numerical simulations of the recurrent network model in 
% the impossible trials. Responses of each network were divided into 5 
% quantiles (quintiles) according to the reaction time (RTs). Within each 
% quintile, the fraction of choices congruent with the overall bias of the 
% network, pBias, was computed and averaged over the different networks 
% (n=200 networks, 500 trials per network). This analysis was performed 
% independently for different values of g. 
% Fig. 7b: CBF for 'IC+drift bias' DDM (fit-based) simulations of the 
% recurrent network simulations.
% Fig. 7c: CBF for the vertical bisection and motor task.



%% Reading and analyzing PPC data for all of g's datasets (8) and the drift+IC model:

PPC_link = 'hddm related code and data\';
dataName_vect = {'g27', 'g28', 'g29', 'g3', ...
    'g31', 'g32', 'g33', 'g34'}; 
dataNamePrint_vect = {'g = 2.7', 'g = 2.8', 'g = 2.9', 'g = 3.0', 'g = 3.1', ...
    'g = 3.2', 'g = 3.3', 'g = 3.4'};
modelName = 'both';
reps = 20; % #samples in PPC per participant and trial:
nObservers = 200; % #networks in data:
nTrials = 500;

% running over g datasets (8):
for dataNum = 1:8
    dataName = dataName_vect(dataNum);
    dataLink = PPC_link + string(dataName) + ...
        '\informative priors 40K\both\ppc_data.csv';
    PPC_data =  csvread( dataLink ,1 ,1 );
    % computing observed and expected probabilities, overall and in 5 
    % RT quantiles:
    [pE_pO_qE1to5_qO1to5_A, mRTbyQuantile_O_A, mRTbyQuantile_E_A ] = ...
        PPC_in_model(PPC_data,reps,nObservers);
    PPC_analysis_struct.( dataName_vect{dataNum} ).both ...
        = pE_pO_qE1to5_qO1to5_A;
    quantileRTm.( dataName_vect{dataNum} ).both.observed ...
        = mRTbyQuantile_O_A;
    quantileRTm.( dataName_vect{dataNum} ).both.expected ...
        = mRTbyQuantile_E_A;
    % RT matrices for observed and PPC-simulated data:
    [RT_observed_A,RT_expected_A,RT_observed_real_A,RT_expected_real_A] = ...
        DT_data(PPC_data,reps,nObservers,nTrials); 
    PPC_RT.( dataName_vect{dataNum} ).both =  RT_expected_A;
    PPC_RT_real.( dataName_vect{dataNum} ).both =  RT_expected_real_A;
    PPC_RT.( dataName_vect{dataNum} ).observed = RT_observed_A;
    PPC_RT_real.( dataName_vect{dataNum} ).observed = RT_observed_real_A;
    PPC_RT_both.( dataName_vect{dataNum} ).RT_observed = RT_observed_A;
    PPC_RT_both.( dataName_vect{dataNum} ).RT_expected = RT_expected_A;
end



%% Computing conditional bias functions:

forModelsCols = ( 1/255 ) * [72 171 72; 234 175 32; 183 82 158; ...
    75 75 75; 82 189 236; 153 153 153; 10 115 184; 198 87 48];

forModelsPlot = {'+-', 'd-', 'o-', 'p-', '*-', 'v-', 's-', '^-'}; 

forNormalizedSlope = linspace( 0.5, 1, 5 );

mean_sem_slope_E = nan(8,2);
mean_sem_slope_O = nan(8,2);

p_max_O_allQs = nan(200,8);
p_max_E_allQs = nan(200,8);

observed_fig = figure;
expected_fig = figure;

for dataNum = 1:8
    
    dataName = dataName_vect(dataNum);

    pE_pO_qE1to5_qO1to5 = ...
        PPC_analysis_struct.( dataName_vect{dataNum} ).both;
   quantileRTmO = quantileRTm.( dataName_vect{dataNum} ).both.observed;
   quantileRTmE = quantileRTm.( dataName_vect{dataNum} ).both.expected;

    % Fig. 7a: CBF for Numerical simulations:
    
    % load relevant data:
    p_max_O_allQs(:,dataNum) = [...
        pE_pO_qE1to5_qO1to5( pE_pO_qE1to5_qO1to5(:,2)>=0.5, 2 ); ...
        1-pE_pO_qE1to5_qO1to5( pE_pO_qE1to5_qO1to5(:,2)<0.5, 2 )];
    p_max_E_allQs(:,dataNum) = [...
        pE_pO_qE1to5_qO1to5( pE_pO_qE1to5_qO1to5(:,1)>=0.5, 1 ); ...
        1-pE_pO_qE1to5_qO1to5( pE_pO_qE1to5_qO1to5(:,1)<0.5, 1 )];
    p_max_observed_up0 = pE_pO_qE1to5_qO1to5( ...
        pE_pO_qE1to5_qO1to5(:,2)>=0.5, 8:12 );
    p_max_observed_down = pE_pO_qE1to5_qO1to5( ...
        pE_pO_qE1to5_qO1to5(:,2)<0.5, 8:12 );
    p_max_observed = [p_max_observed_up0; 1-p_max_observed_down];
    realRT_observed = PPC_RT_real.( dataName_vect{dataNum} ).observed;
    
    % Omitting cases for which the 5th quantile of observed has zero trials. 
    % For the rest, there are exactly 5 trials in a the 1-4 quantiles and 1-5
    % in the 5th.
    p_max_observed = p_max_observed( isnan(p_max_observed(:,5))==0, : );
    
    % plot conditional: 
    figure(observed_fig);
    errorbar( 1:5, mean( p_max_observed ), ...
        std( p_max_observed ) ./ sqrt( size( p_max_observed, 1 ) ), ...
        forModelsPlot{dataNum}, 'Color', forModelsCols(dataNum,:), ...
        'lineWidth', 1); hold on;
    text( 5.2, mean(p_max_observed(:,5)), ...
        dataNamePrint_vect{dataNum}, 'Color', forModelsCols(dataNum,:) );
    v = pca( [ reshape( repmat( forNormalizedSlope, [nObservers,1] ), [], 1 ), ...
        reshape( p_max_observed, [], 1 ) ] );
    bySubSlope = nan(nObservers,1);
    
    % by-network slope:
    for ii = 1:nObservers 
        RT_sub = abs(realRT_observed(ii,:))';
        dec_sub = ( realRT_observed(ii,:) > 0 )';
        p = mean( dec_sub );
        if p < 0.5
            dec_sub = abs(dec_sub - 1);
        end
        RT_p_sub = sortrows([RT_sub, dec_sub]);
        v2 = pca( [(1:length(RT_sub))' RT_p_sub(:,2)] );
        bySubSlope(ii) = v2(2,1)/v2(1,1);
    end
    mean_sem_slope_O(dataNum,:) = [mean(bySubSlope), ...
        std(bySubSlope)/sqrt(nObservers)];
    
    % inset:
    if dataNum == 8
        title( 'Numerical simulations' );
        ylabel('mean p_{max}'); 
        xticks(1:5)
        xticklabels({'0.1','0.3','0.5','0.7','0.9'})
        xlabel('RT (quantiles)'); xlim([1,5.5])
        ylim([0.5,0.9]); pbaspect([1 .75 1])
        box off; 
        ax2 = axes('Position',[0.25 0.2 0.4 0.2]);
        errorbar( 2.7:0.1:3.4, 5*mean_sem_slope_O(:,1), 5*mean_sem_slope_O(:,2), ...
            'Color', 'k', 'Marker', '.', 'lineWidth', 2 );
        xlabel('g'); ylabel('dp/d%'); box off; 
        xlim([2.7,3.4]); xticks(2.7:0.1:3.4); ylim([-2e-3,1e-4]);
        hold on;
    end
    
    
    % Fig. 7b: CBF for 'IC+drift bias' DDM (fit-based) simulations:
    
    % load relevant data:
    p_max_expected_up0 = pE_pO_qE1to5_qO1to5( ...
        pE_pO_qE1to5_qO1to5(:,1)>=0.5, 3:7 );
    p_max_expected_down = pE_pO_qE1to5_qO1to5( ...
        pE_pO_qE1to5_qO1to5(:,1)<0.5, 3:7 );
    p_max_expected = [p_max_expected_up0; 1-p_max_expected_down];
    realRT_expected = PPC_RT_real.( dataName_vect{dataNum} ).both;
        
    % Omitting case with #trials in 5th quantile = 0: 
    p_max_expected = p_max_expected( isnan(p_max_observed(:,5))==0, : );
    
    % plot conditional: 
    figure(expected_fig);
    errorbar( (1:5), mean( p_max_expected ), ...
        std( p_max_expected ) ./ sqrt( size( p_max_expected, 1 ) ), ...
        forModelsPlot{dataNum}, 'Color', forModelsCols(dataNum,:), ...
        'lineWidth', 1); hold on;
    text( 5.2, mean(p_max_expected(:,5)), ...
        dataNamePrint_vect{dataNum}, 'Color', forModelsCols(dataNum,:) );
    v = pca( [ reshape( repmat( forNormalizedSlope, [nObservers,1] ), [], 1 ), ...
        reshape( p_max_expected, [], 1 ) ] );
    bySubSlope = nan(nObservers,1);
    
    % by-network slope:
    for ii = 1:nObservers 
        RT_sub = abs(realRT_expected(ii,:))';
        dec_sub = ( realRT_expected(ii,:) > 0 )';
        p = mean( dec_sub );
        if p < 0.5
            dec_sub = abs(dec_sub - 1);
        end
        RT_p_sub = sortrows([RT_sub, dec_sub]);
        v2 = pca( [(1:length(RT_sub))' RT_p_sub(:,2)] );
        bySubSlope(ii) = v2(2,1)/v2(1,1);
    end
    mean_sem_slope_E(dataNum,:) = [mean(bySubSlope), ...
        std(bySubSlope)/sqrt(nObservers)];
    
    % inset:
    if dataNum == 8
        title( 'DDM-based simulations' );
        ylabel('mean p_{max}'); 
        xticks(1:5)
        xticklabels({'0.1','0.3','0.5','0.7','0.9'})
        xlabel('RT (quantiles)'); xlim([1,5.5])
        ylim([0.5,0.9]); pbaspect([1 .75 1])
        box off; 
        % Slope:
        ax3 = axes('Position',[0.25 0.2 0.4 0.2]);
        errorbar( 2.7:0.1:3.4, 100*mean_sem_slope_E(:,1), 100*mean_sem_slope_E(:,2), ...
            'Color', 'k', 'Marker', '.', 'lineWidth', 2 );
        xlabel('g'); ylabel('dp/d%'); box off; 
        xlim([2.7,3.4]); xticks(2.7:0.1:3.4); ylim([-2e-3,1e-4]); 
        hold on;
    end
    
end 

% color inset datapoints:
gVect = 2.7:0.1:3.4;
for dataNum = 1:8
    % observed:
    figure(observed_fig);
    plot( gVect(dataNum), 5*mean_sem_slope_O(dataNum,1), 'Marker', 'o', ...
        'MarkerFaceColor', forModelsCols(dataNum,:), ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 4 );
    hold on;
    % expected:
    figure(expected_fig);
    plot( gVect(dataNum), 100*mean_sem_slope_E(dataNum,1), 'Marker', 'o', ...
        'MarkerFaceColor', forModelsCols(dataNum,:), ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 4 );
    hold on;
end




%% Fig. 7c: CBF for the vertical bisection and motor task.


% load PPC data:

isPPC = exist('PPC_loop_data.mat'); %
if isPPC == 2
    load('PPC_loop_data.mat') 
else 
    PPC_loop
end

dataName_vect = {'BL', 'g27', 'g3', 'g34', ...
    'motor1', 'motor2', 'motor3', 'motor4', 'motor5', ...
    'motor6', 'motor7', 'motor8', 'motor9', 'motor10'};
modelName_vect = {'none', 'drift', 'ic', 'both'};
modelNum = 4;


figure;


for dataNum = 1:5
    dataName = dataName_vect(dataNum);
    
    % consider only observed and 'IC+drift' DDM posterior simulated data:
    if modelNum == 4 && ( abs(dataNum-3) == 2 )
        
        % load relevant data:
        condCol = [0 0 0];
        modelName = modelName_vect(modelNum);
        pE_pO_qE1to5_qO1to5 = [];
        if dataNum == 5 % if motor then merge all datasets of 10 pairs
            condCol = [.5 .5 .5];
            for dataNumMotor = 5:14
                motorpE_pO_adding = ...
                PPC_analysis_struct.( dataName_vect{dataNumMotor} ).( modelName_vect{modelNum} );
                pE_pO_qE1to5_qO1to5 = [pE_pO_qE1to5_qO1to5; motorpE_pO_adding];
            end
        else
            pE_pO_qE1to5_qO1to5 = ...
            PPC_analysis_struct.( dataName_vect{dataNum} ).( modelName_vect{modelNum} );
        end
        
        % Observed data - numerical simulations of the recurrent net model:
        p_max_observed_up0 = pE_pO_qE1to5_qO1to5( ...
        pE_pO_qE1to5_qO1to5(:,2)>=0.5, 8:12 );
        p_max_observed_down = pE_pO_qE1to5_qO1to5( ...
        pE_pO_qE1to5_qO1to5(:,2)<0.5, 8:12 );
        p_max_observed = [p_max_observed_up0; 1-p_max_observed_down];
        
        % Omitting cases for which the 5th quantile of observed has zero trials.
        % For the rest, there are exactly 5 trials in a the 1-4 quantiles and 1-5
        % in the 5th.
        p_max_observed = p_max_observed( isnan(p_max_observed(:,5))==0, : );
        
        % plot conditional bia function:
        errorbar( 1:5, mean( p_max_observed ), ...
            std( p_max_observed ) ./ sqrt( size( p_max_observed, 1 ) ), ...
            'Color', condCol, 'lineStyle', '-' , 'lineWidth', 1, 'Marker', 'o', ...
            'MarkerFaceColor', condCol, 'MarkerEdgeColor', 'none' );
        hold on;
        
        % 'IC+drift bias' DDM (fit-based) simulations of the recurrent
        % data:
        p_max_expected_up0 = pE_pO_qE1to5_qO1to5( ...
        pE_pO_qE1to5_qO1to5(:,1)>=0.5, 3:7 );
        p_max_expected_down = pE_pO_qE1to5_qO1to5( ...
        pE_pO_qE1to5_qO1to5(:,1)<0.5, 3:7 );
        p_max_expected = [p_max_expected_up0; 1-p_max_expected_down];
        
        % Omitting case with #trials in 5th quantile = 0:
        p_max_expected = p_max_expected( isnan(p_max_observed(:,5))==0, : );
        
        % plot conditional bia function:
        errorbar( (1:5), mean( p_max_expected ), ...
            std( p_max_expected ) ./ sqrt( size( p_max_expected, 1 ) ), ...
            'Color', condCol, 'lineStyle', '--', 'lineWidth', 1, 'Marker', 's', ...
            'MarkerFaceColor', condCol, 'MarkerEdgeColor', 'none' );
        hold on;
    end
    
    ylabel('p_{bias}');
    ggg = gca;
    ggg.YMinorTick = 'on';
    xticks(1:5)
    xticklabels({'QU_1','QU_2','QU_3','QU_4','QU_5'})
    xlabel('RT'); xlim([1,5])
    ylim([0.5,0.9]);
    box off;
    
end
